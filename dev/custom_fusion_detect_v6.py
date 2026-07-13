#!/usr/bin/env python3
"""
Whitelist-based fusion support scanner for BAM files.

What it does
------------
- Reads a TSV whitelist of allowed chromosome-pair templates.
- Supports either fixed-bin targets or midpoint-defined targets.
- Uses a per-row window size when provided; otherwise defaults to the global window size.
- Expands each exact window by padding around the window boundaries.
- Supports padding measured in window units (default) or base pairs.
- Uses targeted BAM fetches for exact and broad windows.
- Calls samtools flagstat for sample-level totals.
- Reports:
    * total_primary_mapped_reads
    * total_chimeric_reads (samtools supplementary alignment count)
    * exact_left_reads / exact_right_reads / exact_support_reads / exact_split_reads
    * broad_left_reads / broad_right_reads / broad_support_reads / broad_split_reads
    * supporting read names and split read names

Important
---------
- Contig names are resolved automatically:
    chr14 <-> 14
- Fusion direction is irrelevant:
  left/right are only used to define the pair and padding.
- Secondary alignments are ignored by default.
- Supplementary alignments are included as ordinary alignments for support counting,
  and reads with an SA tag are also treated as split-supporting.
- The BAM must be indexed because the script uses targeted fetches.
- samtools must be available in PATH.

Whitelist TSV format
--------------------
Required columns:
    cluster, chr_left, chr_right

Target specification:
    Use ONE of the following per row:

    Bin-based (default):
        bin_left, bin_right

    Midpoint-based:
        mid_left_bp, mid_right_bp

Optional columns:
    window_size, pad_left, pad_right

Meaning:
- If mid_left_bp / mid_right_bp are present for a row, they take precedence and bin_* is ignored.
- Otherwise bin_left / bin_right are interpreted as 1-based fixed windows.
- window_size is in base pairs and defaults to the global --window-size (1 Mb by default).
- pad_left and pad_right are interpreted according to --padding-units:
    * windows (default): pad values are counts of windows
    * bp: pad values are base pairs

Examples
--------
Bin-based:
    cluster  chr_left  bin_left  chr_right  bin_right  window_size  pad_left  pad_right
    t4;14    4         22        14         63         1000000      1         1

Midpoint-based:
    cluster  chr_left  mid_left_bp  chr_right  mid_right_bp  window_size  pad_left  pad_right
    t4;14    4         22000000     14        63000000      1000000      250000    250000
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import shutil
import statistics
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pysam


DEFAULT_WINDOW_SIZE = 1_000_000


@dataclass
class Target:
    cluster: str
    chr_left: str
    chr_right: str
    window_size: int
    pad_left: float
    pad_right: float
    padding_units: str
    target_mode: str
    line_no: int

    # Bin or midpoint inputs (stored for reporting)
    bin_left: Optional[int] = None
    bin_right: Optional[int] = None
    mid_left_bp: Optional[int] = None
    mid_right_bp: Optional[int] = None

    exact_left_start: int = 0
    exact_left_end: int = 0
    exact_right_start: int = 0
    exact_right_end: int = 0
    broad_left_start: int = 0
    broad_left_end: int = 0
    broad_right_start: int = 0
    broad_right_end: int = 0


def bin_to_interval(bin_id: int, window_size: int) -> Tuple[int, int]:
    if bin_id < 1:
        raise ValueError(f"bin_id must be >= 1, got {bin_id}")
    start = (bin_id - 1) * window_size + 1
    end = bin_id * window_size
    return start, end


def midpoint_to_interval(mid_bp: int, window_size: int) -> Tuple[int, int]:
    if mid_bp < 1:
        raise ValueError(f"midpoint must be >= 1, got {mid_bp}")
    half = window_size // 2
    start = max(1, mid_bp - half)
    end = start + window_size - 1
    return start, end


def expand_interval(
    start: int,
    end: int,
    pad_value: float,
    window_size: int,
    padding_units: str,
) -> Tuple[int, int]:
    if pad_value < 0:
        raise ValueError(f"pad_value must be >= 0, got {pad_value}")

    if padding_units == "windows":
        pad_bp = pad_value * window_size
    elif padding_units == "bp":
        pad_bp = pad_value
    else:
        raise ValueError(f"Unknown padding_units: {padding_units}")

    pad_bp_int = int(round(pad_bp))
    start = max(1, start - pad_bp_int)
    end = end + pad_bp_int
    return start, end


def build_contig_resolver(bam: pysam.AlignmentFile):
    refs = set(bam.references)
    alias: Dict[str, str] = {}

    for r in refs:
        alias[r] = r
        if r.startswith("chr"):
            alias[r[3:]] = r
        else:
            alias["chr" + r] = r

    def resolve(name: str) -> Optional[str]:
        name = str(name).strip()
        if name in alias:
            return alias[name]
        if name.startswith("chr") and name[3:] in alias:
            return alias[name[3:]]
        if not name.startswith("chr") and ("chr" + name) in alias:
            return alias["chr" + name]
        return None

    return resolve


def parse_optional_int(value: Optional[str], *, default: Optional[int] = None) -> Optional[int]:
    if value is None or value == "":
        return default
    return int(value)


def parse_optional_float(value: Optional[str], *, default: float = 1.0) -> float:
    if value is None or value == "":
        return default
    return float(value)


def load_whitelist(path: str, default_window_size: int) -> List[dict]:
    rows: List[dict] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"cluster", "chr_left", "chr_right"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Whitelist missing columns: " + ", ".join(sorted(missing)))

        for line_no, row in enumerate(reader, start=2):
            rows.append(
                {
                    "cluster": row["cluster"].strip(),
                    "chr_left": row["chr_left"].strip(),
                    "chr_right": row["chr_right"].strip(),
                    "bin_left": parse_optional_int(row.get("bin_left")),
                    "bin_right": parse_optional_int(row.get("bin_right")),
                    "mid_left_bp": parse_optional_int(row.get("mid_left_bp")),
                    "mid_right_bp": parse_optional_int(row.get("mid_right_bp")),
                    "window_size": parse_optional_int(row.get("window_size"), default=default_window_size),
                    "pad_left": parse_optional_float(row.get("pad_left"), default=1.0),
                    "pad_right": parse_optional_float(row.get("pad_right"), default=1.0),
                    "line_no": line_no,
                }
            )
    return rows


def prepare_targets(
    raw_rows: List[dict],
    bam: pysam.AlignmentFile,
    default_padding_units: str,
    logger: logging.Logger,
) -> List[Target]:
    resolve = build_contig_resolver(bam)
    targets: List[Target] = []

    for row in raw_rows:
        chr_left = resolve(row["chr_left"])
        chr_right = resolve(row["chr_right"])

        if chr_left is None:
            raise ValueError(
                f"Could not resolve chr_left='{row['chr_left']}' for whitelist line {row['line_no']}"
            )
        if chr_right is None:
            raise ValueError(
                f"Could not resolve chr_right='{row['chr_right']}' for whitelist line {row['line_no']}"
            )

        window_size = row["window_size"]
        if window_size is None or window_size < 1:
            raise ValueError(
                f"Invalid window_size at whitelist line {row['line_no']}: {window_size}"
            )

        has_midpoints = row["mid_left_bp"] is not None or row["mid_right_bp"] is not None
        has_bins = row["bin_left"] is not None or row["bin_right"] is not None

        if has_midpoints:
            if row["mid_left_bp"] is None or row["mid_right_bp"] is None:
                raise ValueError(
                    f"If either midpoint column is used, both mid_left_bp and mid_right_bp must be present "
                    f"(whitelist line {row['line_no']})"
                )
            exact_left_start, exact_left_end = midpoint_to_interval(row["mid_left_bp"], window_size)
            exact_right_start, exact_right_end = midpoint_to_interval(row["mid_right_bp"], window_size)
            target_mode = "midpoint"
            bin_left = None
            bin_right = None
            mid_left_bp = row["mid_left_bp"]
            mid_right_bp = row["mid_right_bp"]
        elif has_bins:
            if row["bin_left"] is None or row["bin_right"] is None:
                raise ValueError(
                    f"If either bin column is used, both bin_left and bin_right must be present "
                    f"(whitelist line {row['line_no']})"
                )
            exact_left_start, exact_left_end = bin_to_interval(row["bin_left"], window_size)
            exact_right_start, exact_right_end = bin_to_interval(row["bin_right"], window_size)
            target_mode = "bin"
            bin_left = row["bin_left"]
            bin_right = row["bin_right"]
            mid_left_bp = None
            mid_right_bp = None
        else:
            raise ValueError(
                f"Whitelist line {row['line_no']} must provide either bin_left/bin_right or "
                f"mid_left_bp/mid_right_bp"
            )

        broad_left_start, broad_left_end = expand_interval(
            exact_left_start,
            exact_left_end,
            row["pad_left"],
            window_size,
            default_padding_units,
        )
        broad_right_start, broad_right_end = expand_interval(
            exact_right_start,
            exact_right_end,
            row["pad_right"],
            window_size,
            default_padding_units,
        )

        targets.append(
            Target(
                cluster=row["cluster"],
                chr_left=chr_left,
                chr_right=chr_right,
                window_size=window_size,
                pad_left=row["pad_left"],
                pad_right=row["pad_right"],
                padding_units=default_padding_units,
                target_mode=target_mode,
                line_no=row["line_no"],
                bin_left=bin_left,
                bin_right=bin_right,
                mid_left_bp=mid_left_bp,
                mid_right_bp=mid_right_bp,
                exact_left_start=exact_left_start,
                exact_left_end=exact_left_end,
                exact_right_start=exact_right_start,
                exact_right_end=exact_right_end,
                broad_left_start=broad_left_start,
                broad_left_end=broad_left_end,
                broad_right_start=broad_right_start,
                broad_right_end=broad_right_end,
            )
        )

    logger.info("Loaded %d whitelist targets", len(targets))
    return targets


def configure_logging(log_path: Path, verbose: bool) -> logging.Logger:
    logger = logging.getLogger("custom_fusion")
    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    for h in list(logger.handlers):
        logger.removeHandler(h)

    log_path.parent.mkdir(parents=True, exist_ok=True)

    file_handler = logging.FileHandler(log_path, mode="w")
    file_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    file_handler.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s | %(levelname)s | %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    )
    logger.addHandler(file_handler)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO if verbose else logging.ERROR)
    console_handler.setFormatter(logging.Formatter("%(levelname)s | %(message)s"))
    logger.addHandler(console_handler)

    return logger


def derive_default_prefix(bam_path: str) -> str:
    return Path(bam_path).stem


def is_split_alignment(aln) -> bool:
    return aln.is_supplementary or aln.has_tag("SA")


def check_samtools_available() -> None:
    if shutil.which("samtools") is None:
        raise SystemExit(
            "samtools was not found in PATH. Install samtools or add it to PATH before running this script."
        )


def run_samtools_flagstat(bam_path: str, logger: logging.Logger) -> Dict[str, int]:
    check_samtools_available()

    logger.info("Running samtools flagstat")
    try:
        proc = subprocess.run(
            ["samtools", "flagstat", bam_path],
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(
            f"samtools flagstat failed for {bam_path}: {e.stderr or e.stdout or e}"
        ) from e

    total = 0
    primary_mapped = 0
    mapped = 0
    secondary = 0
    supplementary = 0

    line_re = re.compile(r"^\s*(\d+)\s+\+\s+\d+\s+(.+?)(?:\s+\(|$)")
    for raw_line in proc.stdout.splitlines():
        m = line_re.match(raw_line)
        if not m:
            continue
        count = int(m.group(1))
        label = m.group(2).strip()

        if label == "in total":
            total = count
        elif label == "primary mapped":
            primary_mapped = count
        elif label == "mapped":
            mapped = count
        elif label == "secondary":
            secondary = count
        elif label == "supplementary":
            supplementary = count

    if primary_mapped == 0 and mapped:
        primary_mapped = max(mapped - secondary - supplementary, 0)

    logger.info(
        "samtools flagstat summary: total=%d primary_mapped=%d supplementary=%d",
        total,
        primary_mapped,
        supplementary,
    )

    return {
        "total_primary_mapped_reads": primary_mapped,
        "total_chimeric_reads": supplementary,
        "samtools_total_reads": total,
        "samtools_mapped_reads": mapped,
        "samtools_secondary_alignments": secondary,
        "samtools_supplementary_alignments": supplementary,
    }


def collect_reads_in_interval(
    bam: pysam.AlignmentFile,
    chrom: str,
    start_1based: int,
    end_1based: int,
    min_mapq: int,
    keep_secondary: bool,
) -> Tuple[Set[str], Set[str], List[int]]:
    read_names: Set[str] = set()
    supplementary_read_names: Set[str] = set()
    mapqs: List[int] = []

    fetch_start0 = max(0, start_1based - 1)
    fetch_end0 = end_1based

    for aln in bam.fetch(chrom, fetch_start0, fetch_end0):
        if aln.is_unmapped:
            continue
        if not keep_secondary and aln.is_secondary:
            continue
        if aln.mapping_quality < min_mapq:
            continue

        read_names.add(aln.query_name)
        mapqs.append(aln.mapping_quality)

        if is_split_alignment(aln):
            supplementary_read_names.add(aln.query_name)

    return read_names, supplementary_read_names, mapqs


def scan_targets(
    bam: pysam.AlignmentFile,
    targets: List[Target],
    min_mapq: int,
    keep_secondary: bool,
    logger: logging.Logger,
    sample_stats: Dict[str, int],
) -> List[dict]:
    results: List[dict] = []

    logger.info("Scanning targeted windows")
    for t in targets:
        exact_left_reads, exact_left_split_reads, exact_left_mapqs = collect_reads_in_interval(
            bam=bam,
            chrom=t.chr_left,
            start_1based=t.exact_left_start,
            end_1based=t.exact_left_end,
            min_mapq=min_mapq,
            keep_secondary=keep_secondary,
        )

        exact_right_reads, exact_right_split_reads, exact_right_mapqs = collect_reads_in_interval(
            bam=bam,
            chrom=t.chr_right,
            start_1based=t.exact_right_start,
            end_1based=t.exact_right_end,
            min_mapq=min_mapq,
            keep_secondary=keep_secondary,
        )

        exact_supporting = sorted(exact_left_reads & exact_right_reads)
        exact_support_set = set(exact_supporting)
        exact_split_supporting = sorted(
            exact_support_set & (exact_left_split_reads | exact_right_split_reads)
        )
        exact_support_count = len(exact_supporting)
        exact_split_count = len(exact_split_supporting)
        exact_all_mapqs = exact_left_mapqs + exact_right_mapqs
        exact_median_mapq = statistics.median(exact_all_mapqs) if exact_all_mapqs else ""

        broad_left_reads, broad_left_split_reads, broad_left_mapqs = collect_reads_in_interval(
            bam=bam,
            chrom=t.chr_left,
            start_1based=t.broad_left_start,
            end_1based=t.broad_left_end,
            min_mapq=min_mapq,
            keep_secondary=keep_secondary,
        )

        broad_right_reads, broad_right_split_reads, broad_right_mapqs = collect_reads_in_interval(
            bam=bam,
            chrom=t.chr_right,
            start_1based=t.broad_right_start,
            end_1based=t.broad_right_end,
            min_mapq=min_mapq,
            keep_secondary=keep_secondary,
        )

        broad_supporting = sorted(broad_left_reads & broad_right_reads)
        broad_support_set = set(broad_supporting)
        broad_split_supporting = sorted(
            broad_support_set & (broad_left_split_reads | broad_right_split_reads)
        )
        broad_support_count = len(broad_supporting)
        broad_split_count = len(broad_split_supporting)
        broad_all_mapqs = broad_left_mapqs + broad_right_mapqs
        broad_median_mapq = statistics.median(broad_all_mapqs) if broad_all_mapqs else ""

        results.append(
            {
                "cluster": t.cluster,
                "chr_left": t.chr_left,
                "chr_right": t.chr_right,
                "target_mode": t.target_mode,
                "line_no": t.line_no,
                "bin_left": t.bin_left if t.bin_left is not None else "",
                "bin_right": t.bin_right if t.bin_right is not None else "",
                "mid_left_bp": t.mid_left_bp if t.mid_left_bp is not None else "",
                "mid_right_bp": t.mid_right_bp if t.mid_right_bp is not None else "",
                "window_size": t.window_size,
                "pad_left": t.pad_left,
                "pad_right": t.pad_right,
                "padding_units": t.padding_units,
                "total_primary_mapped_reads": sample_stats["total_primary_mapped_reads"],
                "total_chimeric_reads": sample_stats["total_chimeric_reads"],
                "samtools_total_reads": sample_stats["samtools_total_reads"],
                "samtools_mapped_reads": sample_stats["samtools_mapped_reads"],
                "samtools_secondary_alignments": sample_stats["samtools_secondary_alignments"],
                "samtools_supplementary_alignments": sample_stats["samtools_supplementary_alignments"],
                "exact_left_interval": f"{t.chr_left}:{t.exact_left_start}-{t.exact_left_end}",
                "exact_right_interval": f"{t.chr_right}:{t.exact_right_start}-{t.exact_right_end}",
                "broad_left_interval": f"{t.chr_left}:{t.broad_left_start}-{t.broad_left_end}",
                "broad_right_interval": f"{t.chr_right}:{t.broad_right_start}-{t.broad_right_end}",
                "exact_left_reads": len(exact_left_reads),
                "exact_right_reads": len(exact_right_reads),
                "exact_support_reads": exact_support_count,
                "exact_split_reads": exact_split_count,
                "exact_median_mapq": exact_median_mapq,
                "exact_supporting_read_names": ";".join(exact_supporting),
                "exact_split_read_names": ";".join(exact_split_supporting),
                "exact_left_supplementary_reads": len(exact_left_split_reads),
                "exact_right_supplementary_reads": len(exact_right_split_reads),
                "broad_left_reads": len(broad_left_reads),
                "broad_right_reads": len(broad_right_reads),
                "broad_support_reads": broad_support_count,
                "broad_split_reads": broad_split_count,
                "broad_median_mapq": broad_median_mapq,
                "broad_supporting_read_names": ";".join(broad_supporting),
                "broad_split_read_names": ";".join(broad_split_supporting),
                "broad_left_supplementary_reads": len(broad_left_split_reads),
                "broad_right_supplementary_reads": len(broad_right_split_reads),
            }
        )

        logger.debug(
            "%s | mode=%s exact_support=%d exact_split=%d broad_support=%d broad_split=%d",
            t.cluster,
            t.target_mode,
            exact_support_count,
            exact_split_count,
            broad_support_count,
            broad_split_count,
        )

    return results


def make_output_paths(outdir: Path, out_prefix: str) -> Tuple[Path, Path]:
    log_path = outdir / f"{out_prefix}.log"
    tsv_path = outdir / f"{out_prefix}.fusion_whitelist.tsv"
    return log_path, tsv_path


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Whitelist-based fusion support scanner from BAM using fixed or midpoint-defined genomic windows. "
            "Fusion direction is treated as irrelevant: the script only checks whether a read hits both sides "
            "of a user-specified chromosome-pair template."
        )
    )
    parser.add_argument("--bam", required=True, help="Input BAM file.")
    parser.add_argument(
        "--whitelist",
        required=True,
        help=(
            "TSV whitelist with required columns: cluster, chr_left, chr_right. "
            "Use either bin_left/bin_right or mid_left_bp/mid_right_bp per row. "
            "Optional columns: window_size, pad_left, pad_right."
        ),
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help=(
            "Output directory. Default: current working directory. "
            "The directory is created if needed."
        ),
    )
    parser.add_argument(
        "--out-prefix",
        default=None,
        help=(
            "Output filename prefix. Default: BAM basename without the .bam suffix. "
            "Example: path/sample.bam -> sample.log and sample.fusion_whitelist.tsv"
        ),
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=DEFAULT_WINDOW_SIZE,
        help="Default window size in bp for rows that do not provide window_size (default: 1000000).",
    )
    parser.add_argument(
        "--padding-units",
        choices=("windows", "bp"),
        default="windows",
        help=(
            "How to interpret pad_left and pad_right. "
            "windows = padding is measured in window counts (default). "
            "bp = padding is measured directly in base pairs."
        ),
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=20,
        help="Minimum MAPQ to count support (default: 20).",
    )
    parser.add_argument(
        "--keep-secondary",
        action="store_true",
        help="Keep secondary alignments (default: ignore them).",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Write a more detailed process log.",
    )
    args = parser.parse_args()

    if args.out_prefix is None:
        args.out_prefix = derive_default_prefix(args.bam)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    log_path, output_path = make_output_paths(outdir, args.out_prefix)
    logger = configure_logging(log_path, args.verbose)

    logger.info("Starting custom fusion whitelist scan")
    logger.info("Input BAM: %s", args.bam)
    logger.info("Whitelist: %s", args.whitelist)
    logger.info("Output directory: %s", outdir)
    logger.info("Output prefix: %s", args.out_prefix)
    logger.info("Output TSV: %s", output_path)
    logger.info("Log file: %s", log_path)
    logger.info("Default window size: %d bp", args.window_size)
    logger.info("Padding units: %s", args.padding_units)
    logger.info("Minimum MAPQ: %d", args.min_mapq)
    logger.info("Secondary alignments kept: %s", args.keep_secondary)
    logger.info("Fusion direction is treated as irrelevant")

    check_samtools_available()

    bam = pysam.AlignmentFile(args.bam, "rb")
    try:
        if not bam.has_index():
            raise RuntimeError(
                f"BAM index not found for {args.bam}. Create one with samtools index before running this script."
            )

        sample_stats = run_samtools_flagstat(args.bam, logger)
        raw_rows = load_whitelist(args.whitelist, args.window_size)
        targets = prepare_targets(raw_rows, bam, args.padding_units, logger)
        results = scan_targets(
            bam=bam,
            targets=targets,
            min_mapq=args.min_mapq,
            keep_secondary=args.keep_secondary,
            logger=logger,
            sample_stats=sample_stats,
        )

        fieldnames = [
            "cluster",
            "chr_left",
            "chr_right",
            "target_mode",
            "line_no",
            "bin_left",
            "bin_right",
            "mid_left_bp",
            "mid_right_bp",
            "window_size",
            "pad_left",
            "pad_right",
            "padding_units",
            "total_primary_mapped_reads",
            "total_chimeric_reads",
            "samtools_total_reads",
            "samtools_mapped_reads",
            "samtools_secondary_alignments",
            "samtools_supplementary_alignments",
            "exact_left_interval",
            "exact_right_interval",
            "broad_left_interval",
            "broad_right_interval",
            "exact_left_reads",
            "exact_right_reads",
            "exact_support_reads",
            "exact_split_reads",
            "exact_median_mapq",
            "exact_supporting_read_names",
            "exact_split_read_names",
            "exact_left_supplementary_reads",
            "exact_right_supplementary_reads",
            "broad_left_reads",
            "broad_right_reads",
            "broad_support_reads",
            "broad_split_reads",
            "broad_median_mapq",
            "broad_supporting_read_names",
            "broad_split_read_names",
            "broad_left_supplementary_reads",
            "broad_right_supplementary_reads",
        ]

        with open(output_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
            writer.writeheader()
            for row in results:
                writer.writerow(row)

        logger.info("Wrote %d rows to %s", len(results), output_path)
        logger.info("Done")

    except Exception:
        logger.exception("Pipeline failed")
        raise
    finally:
        bam.close()


if __name__ == "__main__":
    main()
