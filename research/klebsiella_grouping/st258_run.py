"""Reproduce the full cached ST258 alignment and bounded-memory distance evidence."""

import argparse
import csv
import json
from pathlib import Path
import resource
import shutil
import time

from chronoclade.cgmlst.locus_alignment import align_sequences, load_alignment, file_hash
from chronoclade.cgmlst.aligned_matrix import build_aligned_matrix
from chronoclade.context_partitions import profile_records
from chronoclade.datasets import load_dataset
from research.klebsiella_grouping.st258_catalogue import digest_object


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset", type=Path, default=Path("tmp/st258-native-e2e/prepared/dataset.json")
    )
    parser.add_argument(
        "--catalogue", type=Path, default=Path("tmp/science-st258/catalogue/catalogue.json")
    )
    parser.add_argument("--output", type=Path, default=Path("tmp/science-st258"))
    parser.add_argument("--phase", choices=["align", "matrix", "all"], default="all")
    args = parser.parse_args()
    start = time.monotonic()
    dataset = load_dataset(args.dataset)
    if len(dataset.profiles) != 1:
        raise ValueError("One preserved profile matrix required")
    records = profile_records(dataset)
    matrix = dataset.profiles[0]
    if tuple(row["sample_id"] for row in records) != matrix.sample_ids:
        by_id = {row["sample_id"]: row for row in records}
        records = [by_id[ident] for ident in matrix.sample_ids]
    catalogue = json.loads(args.catalogue.read_text())
    dna_path = args.catalogue.parent / catalogue["artifact"]["path"]
    if file_hash(dna_path) != catalogue["artifact"]["sha256"]:
        raise ValueError("Vetted allele DNA checksum changed")
    bridge = json.loads((args.catalogue.parent / catalogue["bridge"]["path"]).read_text())
    if (
        digest_object(bridge) != catalogue["bridge"]["canonical_content_sha256"]
        or bridge["original_dataset_id"] != dataset.dataset_id
        or bridge["comparison_loci"] != list(matrix.catalogue.loci)
        or bridge["historical_fingerprints_verified"] != len(records)
    ):
        raise ValueError("Vetted DNA bridge does not match the original dataset")
    with dna_path.open(newline="") as stream:
        values = list(csv.DictReader(stream))
    sequences = {(row["locus"], row["allele"]): row["dna"] for row in values}
    if len(sequences) != len(values):
        raise ValueError("Duplicate vetted allele records")
    print(
        f"Panel: {len(records)} samples, {len(matrix.catalogue.loci)} observed loci, {len(sequences)} resolved alleles",
        flush=True,
    )

    def progress(phase, done, total):
        print(f"{phase}: {done}/{total}; elapsed {time.monotonic() - start:.1f}s", flush=True)

    directory = args.output / "alignment"
    if args.phase in {"align", "all"}:
        if (directory / "alignment.json").is_file():
            raise ValueError("Completed alignment exists; use --phase matrix to verify/reuse")
        executable = shutil.which("mafft")
        if executable is None:
            raise ValueError("Run through the Pixi cgmlst/default environment for MAFFT")
        alignment = align_sequences(
            sequences,
            records,
            directory,
            executable,
            scope=catalogue["scope"],
            provenance={
                "vetted_catalogue": str(args.catalogue.resolve()),
                "vetted_catalogue_sha256": file_hash(args.catalogue),
                "alleles_sha256": file_hash(dna_path),
                "bridge_sha256": file_hash(args.catalogue.parent / catalogue["bridge"]["path"]),
                "original_dataset_id": dataset.dataset_id,
                "original_provider_scheme": matrix.catalogue.scheme_id,
                "original_provider_version": matrix.catalogue.scheme_version,
                "historical_provider_database_fingerprint": matrix.catalogue.database_sha256,
                "observed_locus_catalogue_complete": matrix.catalogue.complete,
                "bridge_scope_is_not_historical_database_identity": True,
            },
            materialize_samples=False,
            export_concatenation=True,
            progress=progress,
        )
    else:
        alignment = load_alignment(directory, records, materialize_samples=False)
    if args.phase in {"matrix", "all"}:
        if (args.output / "nucleotide-matrices" / "matrix.json").is_file():
            raise ValueError("Completed matrices exist; preserve them and choose a new output")
        result = build_aligned_matrix(
            matrix,
            alignment,
            args.output / "nucleotide-matrices",
            min_loci=597,
            block_size=128,
            progress=progress,
        )
        print(f"Matrices complete: {len(result.sample_ids)} samples", flush=True)
    receipt = {
        "phase": args.phase,
        "elapsed_seconds": time.monotonic() - start,
        "peak_main_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "sample_count": len(records),
        "dataset_id": dataset.dataset_id,
        "driver_sha256": file_hash(Path(__file__)),
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / f"{args.phase}-execution.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt), flush=True)


if __name__ == "__main__":
    main()
