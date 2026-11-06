#!/usr/bin/env python3
"""Select a reproducible random OpenDEL sample stratified by qDOS size.

The default target is 26,000 raw SMILES strings. Quotas are calculated
with the largest-remainder method from the original qDOS enumeration counts
supplied for this dataset. Input chunks are expected to have been produced by
deduplicate_opendel_chunks.py.

Only chunk_*.smi files are read. NPY files are never opened.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import heapq
import math
import os
import random
import re
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


# Counts copied from `wc` on OpeDELLibrary/smi/qDOS*.enumeration.smi.
SOURCE_COUNTS = {
    "qDOS11": 77_051_288,
    "qDOS13": 4_287_380,
    "qDOS14": 30_342_972,
    "qDOS18_3": 92_991_800,
    "qDOS19_1": 27_810_720,
    "qDOS19_2": 28_226_000,
    "qDOS23": 77_522_144,
    "qDOS24": 9_424_800,
    "qDOS30": 47_389_944,
    "qDOS42_2": 44_928_000,
    "qDOS4": 234_685_620,
    "qDOS9": 223_650_368,
}

# The molecule IDs expose the base family (for example qDOS9-19-97-1).
# Split source files such as qDOS19_1 and qDOS19_2 are therefore aggregated
# into qDOS19. Their total statistical weight is preserved exactly.
QDOS_ID_RE = re.compile(rb"^qDOS(\d+)")
CHUNK_RE = re.compile(r"^(.*?)(\d+)([^\d]*)$")


def natural_key(path: Path) -> tuple[str, int, str]:
    match = CHUNK_RE.match(path.name)
    if not match:
        return path.name, -1, ""
    return match.group(1), int(match.group(2)), match.group(3)


def family_sort_key(family: str) -> int:
    return int(family.removeprefix("qDOS"))


def aggregate_family_counts() -> dict[str, int]:
    totals: Counter[str] = Counter()
    for source_name, count in SOURCE_COUNTS.items():
        match = re.match(r"qDOS(\d+)", source_name)
        if not match:  # defensive: constants above are controlled
            raise ValueError(f"cannot classify source count {source_name!r}")
        totals[f"qDOS{match.group(1)}"] += count
    return dict(totals)


def calculate_quotas(
    counts: dict[str, int], target: int
) -> tuple[dict[str, int], dict[str, float]]:
    """Allocate integer quotas that sum exactly to target."""
    total = sum(counts.values())
    ideals = {family: target * count / total for family, count in counts.items()}
    quotas = {family: math.floor(ideal) for family, ideal in ideals.items()}
    remaining = target - sum(quotas.values())

    remainder_order = sorted(
        counts,
        key=lambda family: (
            -(ideals[family] - quotas[family]),
            family_sort_key(family),
        ),
    )
    for family in remainder_order[:remaining]:
        quotas[family] += 1

    if sum(quotas.values()) != target:
        raise AssertionError("quota allocation does not sum to target")
    return quotas, ideals


def parse_record(
    line: bytes, source: Path, line_number: int, accepted_families: set[str]
) -> tuple[bytes, bytes, str]:
    record = line.rstrip(b"\r\n")
    if not record:
        raise ValueError(f"{source}:{line_number}: empty line")

    if b"\t" in record:
        smiles, molecule_id = record.split(b"\t", 1)
    else:
        fields = record.split(None, 1)
        if len(fields) != 2:
            raise ValueError(
                f"{source}:{line_number}: expected two columns (SMILES and ID)"
            )
        smiles, molecule_id = fields

    molecule_id = molecule_id.strip()
    if not smiles or not molecule_id:
        raise ValueError(f"{source}:{line_number}: empty SMILES or molecule ID")

    match = QDOS_ID_RE.match(molecule_id)
    if not match:
        shown = molecule_id[:100].decode("utf-8", errors="replace")
        raise ValueError(
            f"{source}:{line_number}: ID does not start with qDOS<number>: {shown!r}"
        )
    family = f"qDOS{int(match.group(1))}"
    if family not in accepted_families:
        shown = molecule_id[:100].decode("utf-8", errors="replace")
        raise ValueError(
            f"{source}:{line_number}: unexpected qDOS family {family}: {shown!r}"
        )
    return smiles, molecule_id, family


@dataclass(frozen=True)
class Record:
    random_key: int
    family: str
    smiles: bytes
    molecule_id: bytes
    chunk: str
    line_number: int


@dataclass(frozen=True)
class WorkerJob:
    worker_number: int
    chunks: tuple[str, ...]
    quotas: dict[str, int]
    seed: int


@dataclass
class WorkerResult:
    worker_number: int
    candidates: dict[str, list[Record]]
    family_rows: dict[str, int]
    input_rows: int
    chunk_rows: dict[str, int]
    elapsed_seconds: float


def _file_seed(seed: int, filename: str) -> int:
    material = f"{seed}\0{filename}".encode("utf-8")
    return int.from_bytes(hashlib.blake2b(material, digest_size=16).digest(), "big")


def _offer(
    heap: list[tuple[int, int, Record]],
    quota: int,
    serial: int,
    record: Record,
) -> None:
    if quota == 0:
        return
    entry = (-record.random_key, serial, record)
    if len(heap) < quota:
        heapq.heappush(heap, entry)
    elif record.random_key < -heap[0][0]:
        heapq.heapreplace(heap, entry)


def scan_partition(job: WorkerJob) -> WorkerResult:
    started = time.monotonic()
    accepted_families = set(job.quotas)
    heaps: dict[str, list[tuple[int, int, Record]]] = {
        family: [] for family in job.quotas
    }
    family_rows: Counter[str] = Counter()
    chunk_rows: dict[str, int] = {}
    input_rows = 0
    serial = 0

    for chunk_number, chunk_string in enumerate(job.chunks, start=1):
        chunk = Path(chunk_string)
        rng = random.Random(_file_seed(job.seed, chunk.name))
        this_chunk_rows = 0

        with chunk.open("rb", buffering=1024 * 1024) as handle:
            for line_number, line in enumerate(handle, start=1):
                input_rows += 1
                this_chunk_rows += 1
                smiles, molecule_id, family = parse_record(
                    line, chunk, line_number, accepted_families
                )

                family_rows[family] += 1
                random_key = rng.getrandbits(128)
                heap = heaps[family]
                quota = job.quotas[family]

                # Construct a Record only if this row enters the reservoir.
                if quota and (len(heap) < quota or random_key < -heap[0][0]):
                    serial += 1
                    record = Record(
                        random_key=random_key,
                        family=family,
                        smiles=smiles,
                        molecule_id=molecule_id,
                        chunk=chunk.name,
                        line_number=line_number,
                    )
                    _offer(heap, quota, serial, record)

        chunk_rows[chunk.name] = this_chunk_rows
        print(
            f"[worker {job.worker_number}] "
            f"{chunk_number}/{len(job.chunks)} {chunk.name}: "
            f"{this_chunk_rows:,} rows",
            flush=True,
        )

    candidates = {
        family: [entry[2] for entry in heap] for family, heap in heaps.items()
    }
    return WorkerResult(
        worker_number=job.worker_number,
        candidates=candidates,
        family_rows=dict(family_rows),
        input_rows=input_rows,
        chunk_rows=chunk_rows,
        elapsed_seconds=time.monotonic() - started,
    )


def merge_worker_results(
    results: list[WorkerResult], quotas: dict[str, int]
) -> tuple[list[Record], Counter[str], int, dict[str, int]]:
    available: Counter[str] = Counter()
    input_rows = 0
    chunk_rows: dict[str, int] = {}
    candidates: dict[str, list[Record]] = defaultdict(list)

    for result in results:
        available.update(result.family_rows)
        input_rows += result.input_rows
        chunk_rows.update(result.chunk_rows)
        for family, records in result.candidates.items():
            candidates[family].extend(records)

    selected: list[Record] = []
    for family in sorted(quotas, key=family_sort_key):
        quota = quotas[family]
        if available[family] < quota:
            raise RuntimeError(
                f"{family}: only {available[family]:,} available rows "
                f"for a quota of {quota:,}"
            )

        family_selected = sorted(
            candidates[family], key=lambda item: item.random_key
        )[:quota]
        if len(family_selected) != quota:
            raise RuntimeError(
                f"{family}: could not fill quota "
                f"({len(family_selected):,}/{quota:,})"
            )
        selected.extend(family_selected)

    return selected, available, input_rows, chunk_rows


def _temporary_path(final_path: Path) -> Path:
    return final_path.with_name(f".{final_path.name}.{os.getpid()}.tmp")


def _prepare_output(path: Path, overwrite: bool) -> Path:
    if path.exists() and not overwrite:
        raise FileExistsError(f"output already exists: {path}; use --overwrite")
    temporary = _temporary_path(path)
    if temporary.exists():
        temporary.unlink()
    return temporary


def write_outputs(
    output_dir: Path,
    selected: list[Record],
    family_counts: dict[str, int],
    ideals: dict[str, float],
    quotas: dict[str, int],
    available: Counter[str],
    input_rows: int,
    chunk_rows: dict[str, int],
    seed: int,
    input_dir: Path,
    overwrite: bool,
) -> tuple[Path, Path, Path, Path]:
    target = len(selected)
    stem = f"opendel_random_{target}"
    smi_path = output_dir / f"{stem}.smi"
    manifest_path = output_dir / f"{stem}_manifest.tsv"
    quotas_path = output_dir / f"{stem}_quotas.tsv"
    summary_path = output_dir / f"{stem}_summary.txt"
    paths = [smi_path, manifest_path, quotas_path, summary_path]
    temporary_paths = [_prepare_output(path, overwrite) for path in paths]

    # Randomize presentation order independently of the bottom-key selection.
    ordered = list(selected)
    random.Random(seed ^ 0x5EED_26_000).shuffle(ordered)

    selected_by_family = Counter(record.family for record in ordered)
    total_source = sum(family_counts.values())
    total_available = sum(available.values())

    try:
        with temporary_paths[0].open("wb", buffering=1024 * 1024) as handle:
            for record in ordered:
                handle.write(record.smiles + b"\t" + record.molecule_id + b"\n")

        with temporary_paths[1].open("wb", buffering=1024 * 1024) as handle:
            handle.write(
                b"sample_index\tqdos_family\tsmiles\tmolecule_id\t"
                b"source_chunk\tsource_line_1based\trandom_key_hex\n"
            )
            for sample_index, record in enumerate(ordered, start=1):
                fields = [
                    str(sample_index).encode(),
                    record.family.encode(),
                    record.smiles,
                    record.molecule_id,
                    record.chunk.encode(),
                    str(record.line_number).encode(),
                    f"{record.random_key:032x}".encode(),
                ]
                handle.write(b"\t".join(fields) + b"\n")

        with temporary_paths[2].open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "qdos_family",
                    "original_count",
                    "original_fraction",
                    "ideal_quota",
                    "integer_quota",
                    "available_deduplicated_chunk_rows",
                    "selected",
                ]
            )
            for family in sorted(family_counts, key=family_sort_key):
                writer.writerow(
                    [
                        family,
                        family_counts[family],
                        f"{family_counts[family] / total_source:.12f}",
                        f"{ideals[family]:.12f}",
                        quotas[family],
                        available[family],
                        selected_by_family[family],
                    ]
                )
            writer.writerow(
                [
                    "TOTAL",
                    total_source,
                    "1.000000000000",
                    f"{sum(ideals.values()):.12f}",
                    sum(quotas.values()),
                    total_available,
                    len(ordered),
                ]
            )

        with temporary_paths[3].open("w", encoding="utf-8") as handle:
            handle.write(f"input_dir\t{input_dir}\n")
            handle.write(f"chunks_read\t{len(chunk_rows)}\n")
            handle.write(f"seed\t{seed}\n")
            handle.write(f"target\t{target}\n")
            handle.write(f"original_qdos_rows\t{total_source}\n")
            handle.write(f"physical_chunk_rows_read\t{input_rows}\n")
            handle.write(f"deduplicated_chunk_rows_seen\t{total_available}\n")
            handle.write(
                f"deduplicated_rows_vs_original_percent\t"
                f"{100 * total_available / total_source:.6f}\n"
            )
            handle.write("qdos_grouping\tbase_qDOS_number\n")

        for temporary, final in zip(temporary_paths, paths):
            os.replace(temporary, final)
    except Exception:
        for temporary in temporary_paths:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass
        raise

    return smi_path, manifest_path, quotas_path, summary_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Randomly select OpenDEL SMILES from deduplicated chunks, with quotas proportional "
            "to the original qDOS enumeration sizes."
        )
    )
    parser.add_argument("--input-dir", type=Path, default=Path("."))
    parser.add_argument("--output-dir", type=Path, default=Path("../OpenDEL_random_26k"))
    parser.add_argument("--pattern", default="chunk_*.smi")
    parser.add_argument("--target", type=int, default=26_000)
    parser.add_argument("--seed", type=int, default=20_260_807)
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="parallel readers; use 1 on the original chunks to limit RAM",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.target < 1:
        raise SystemExit("--target must be at least 1")
    if args.workers < 1:
        raise SystemExit("--workers must be at least 1")

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    if not input_dir.is_dir():
        raise SystemExit(f"input directory does not exist: {input_dir}")
    if input_dir == output_dir:
        raise SystemExit("input and output directories must be different")

    chunks = sorted(input_dir.glob(args.pattern), key=natural_key)
    if not chunks:
        raise SystemExit(f"no files match {args.pattern!r} in {input_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    family_counts = aggregate_family_counts()
    quotas, ideals = calculate_quotas(family_counts, args.target)

    print("qDOS quotas:")
    for family in sorted(quotas, key=family_sort_key):
        print(f"  {family}: {quotas[family]:,}")
    print(f"  TOTAL: {sum(quotas.values()):,}")
    print(f"Reading {len(chunks):,} SMI chunks with {args.workers} worker(s)...")

    worker_count = min(args.workers, len(chunks))
    partitions = [chunks[index::worker_count] for index in range(worker_count)]
    jobs = [
        WorkerJob(
            worker_number=index + 1,
            chunks=tuple(str(path) for path in partition),
            quotas=quotas,
            seed=args.seed,
        )
        for index, partition in enumerate(partitions)
    ]

    if worker_count == 1:
        results = [scan_partition(jobs[0])]
    else:
        results = []
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=worker_count
        ) as executor:
            futures = [executor.submit(scan_partition, job) for job in jobs]
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())

    selected, available, input_rows, chunk_rows = (
        merge_worker_results(results, quotas)
    )
    smi_path, manifest_path, quotas_path, summary_path = write_outputs(
        output_dir=output_dir,
        selected=selected,
        family_counts=family_counts,
        ideals=ideals,
        quotas=quotas,
        available=available,
        input_rows=input_rows,
        chunk_rows=chunk_rows,
        seed=args.seed,
        input_dir=input_dir,
        overwrite=args.overwrite,
    )

    print(f"Done: {len(selected):,} molecules selected.")
    print(f"SMI:      {smi_path}")
    print(f"Manifest: {manifest_path}")
    print(f"Quotas:   {quotas_path}")
    print(f"Summary:  {summary_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except BrokenPipeError:
        try:
            sys.stdout.close()
        finally:
            raise SystemExit(1)

