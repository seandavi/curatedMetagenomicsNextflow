#!/usr/bin/env python
"""Build a per-sample manifest.json compiling provenance + read accounting.

Reads are profiled with a single streaming pass (no external tools), so this
runs inside the base image without adding a dependency. Tool versions are
merged from the per-process versions.yml files staged alongside this script.

See docs/adr/0005-per-sample-manifest.md for the rationale.
"""

import argparse
import glob
import gzip
import json
import re
import sys

VERSION_LINE = re.compile(r"^([A-Za-z0-9_.\-]+):\s*(.+)$")


def _open(path):
    """Open a FASTQ whether or not it is gzip-compressed."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def fastq_stats(path):
    """Single-pass read/base/length/GC statistics over a FASTQ file.

    Returns zeros for an empty or missing file so stub runs and edge cases do
    not break manifest generation.
    """
    reads = 0
    bases = 0
    gc = 0
    min_len = None
    max_len = None
    # Length histogram keeps memory bounded (read lengths are small integers)
    # while still allowing an exact median.
    length_hist = {}

    try:
        with _open(path) as fh:
            for i, line in enumerate(fh):
                if i % 4 != 1:  # only sequence lines
                    continue
                seq = line.strip()
                n = len(seq)
                reads += 1
                bases += n
                gc += seq.count("G") + seq.count("C") + seq.count("g") + seq.count("c")
                length_hist[n] = length_hist.get(n, 0) + 1
                if min_len is None or n < min_len:
                    min_len = n
                if max_len is None or n > max_len:
                    max_len = n
    except FileNotFoundError:
        pass

    median_len = None
    if reads:
        target = reads // 2
        cumulative = 0
        for length in sorted(length_hist):
            cumulative += length_hist[length]
            if cumulative > target:
                median_len = length
                break

    return {
        "number_reads": reads,
        "number_bases": bases,
        "minimum_read_length": min_len if min_len is not None else 0,
        "maximum_read_length": max_len if max_len is not None else 0,
        "median_read_length": median_len if median_len is not None else 0,
        "mean_read_length": round(bases / reads, 2) if reads else 0,
        "gc_percent": round(100.0 * gc / bases, 2) if bases else 0,
    }


def merge_versions(pattern="versions_*.yml"):
    """Merge per-process versions.yml files into a single tool->version map."""
    versions = {}
    for path in sorted(glob.glob(pattern)):
        try:
            with open(path) as fh:
                for line in fh:
                    stripped = line.strip()
                    if not stripped or stripped == "versions:":
                        continue
                    match = VERSION_LINE.match(stripped)
                    if match:
                        versions[match.group(1)] = match.group(2).strip()
        except OSError:
            continue
    return versions


def none_if_blank(value):
    return value if value else None


def as_bool(value):
    return str(value).lower() == "true"


def read_accounting(raw, decontaminated):
    """Raw vs. host-decontaminated stats plus derived survival fractions."""
    block = {
        "raw": raw,
        "decontaminated": decontaminated,
    }
    if raw["number_reads"]:
        block["reads_surviving_fraction"] = round(
            decontaminated["number_reads"] / raw["number_reads"], 4
        )
    if raw["number_bases"]:
        block["bases_surviving_fraction"] = round(
            decontaminated["number_bases"] / raw["number_bases"], 4
        )
    return block


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True)
    p.add_argument("--pipeline-name", default="")
    p.add_argument("--pipeline-version", default="")
    p.add_argument("--nextflow-version", default="")
    p.add_argument("--container", default="")
    p.add_argument("--command-line", default="")
    p.add_argument("--run-name", default="")
    p.add_argument("--session-id", default="")
    p.add_argument("--git-commit", default="")
    p.add_argument("--git-revision", default="")
    p.add_argument("--start-time", default="")
    p.add_argument("--input-mode", default="")
    p.add_argument("--input-ids", default="", help="comma-separated accessions or paths")
    p.add_argument("--metaphlan-index", default="")
    p.add_argument("--store-dir", default="")
    p.add_argument("--skip-rarefied", default="false")
    p.add_argument("--rarefy-reads", type=int, default=0)
    p.add_argument("--rarefy-seed", type=int, default=0)
    p.add_argument("--skip-humann", default="true")
    p.add_argument("--skip-gtdb", default="false")
    p.add_argument("--raw-fastq", required=True)
    p.add_argument("--kneaddata-fastq", required=True)
    p.add_argument("--output", default="manifest.json")
    args = p.parse_args()

    raw = fastq_stats(args.raw_fastq)
    decontaminated = fastq_stats(args.kneaddata_fastq)

    manifest = {
        "sample_id": args.sample,
        "provenance": {
            "pipeline_name": none_if_blank(args.pipeline_name),
            "pipeline_version": none_if_blank(args.pipeline_version),
            "nextflow_version": none_if_blank(args.nextflow_version),
            "container": none_if_blank(args.container),
            "command_line": none_if_blank(args.command_line),
            "run_name": none_if_blank(args.run_name),
            "session_id": none_if_blank(args.session_id),
            "git_commit": none_if_blank(args.git_commit),
            "git_revision": none_if_blank(args.git_revision),
            "start_time": none_if_blank(args.start_time),
            "input_mode": none_if_blank(args.input_mode),
            "input_ids": [x for x in args.input_ids.split(",") if x],
        },
        "parameters": {
            "metaphlan_index": none_if_blank(args.metaphlan_index),
            "store_dir": none_if_blank(args.store_dir),
            "skip_humann": as_bool(args.skip_humann),
            "skip_gtdb": as_bool(args.skip_gtdb),
            "skip_rarefied": as_bool(args.skip_rarefied),
        },
        "read_accounting": read_accounting(raw, decontaminated),
        "software_versions": merge_versions(),
    }

    if not as_bool(args.skip_rarefied):
        manifest["parameters"]["rarefaction"] = {
            "rarefy_reads": args.rarefy_reads,
            "rarefy_seed": args.rarefy_seed,
        }

    with open(args.output, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=False)
        fh.write("\n")

    sys.stderr.write("Wrote manifest for sample {}\n".format(args.sample))


if __name__ == "__main__":
    main()
