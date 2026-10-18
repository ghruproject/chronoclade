"""Align frozen cgMLST allele DNA per locus and reconstruct sample supermatrices.

Missing evidence is N, aligned indels are '-', and substitutions count only
jointly A/C/G/T columns. No positional fallback is available.
"""

from dataclasses import dataclass
from functools import lru_cache
import hashlib
from itertools import combinations
import json
from pathlib import Path
import subprocess

import numpy as np

from chronoclade.cgmlst.dna import read_catalogue
from chronoclade.profile_types import _allele, _profile

MAFFT_OPTIONS = ("--thread", "1", "--inputorder", "--auto", "--nuc", "--quiet")


def file_hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def read_fasta(text):
    """Reject duplicate identifiers, empty records and data before a header."""
    records, ident = {}, None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            ident = line[1:]
            if not ident or any(c.isspace() for c in ident) or ident in records:
                raise ValueError("Invalid or repeated alignment FASTA identifier")
            records[ident] = ""
        elif ident is None:
            raise ValueError("Alignment FASTA sequence precedes its identifier")
        else:
            records[ident] += line.upper()
    if any(not seq for seq in records.values()):
        raise ValueError("Alignment FASTA contains an empty sequence")
    return records


def fasta_text(records):
    return "".join(f">{ident}\n{seq}\n" for ident, seq in records.items())


def validate_msa(aligned, original):
    if not aligned or set(aligned) != set(original):
        raise ValueError("Aligned identifiers differ from input identifiers")
    if len({len(seq) for seq in aligned.values()}) != 1:
        raise ValueError("Aligned sequences have unequal widths")
    for ident, sequence in aligned.items():
        if not sequence or set(sequence) - set("ACGT-"):
            raise ValueError("Alignment contains empty or unsupported DNA")
        if sequence.replace("-", "") != original[ident]:
            raise ValueError("Alignment changed an input ungapped sequence")
    return len(next(iter(aligned.values())))


@lru_cache(maxsize=None)
def _scalar_stats(left, right):
    """Independent scalar arithmetic, with direction-specific indel tracts."""
    snps = sites = gap_base = gap_gap = missing = tracts = previous = 0
    for a, b in zip(left, right, strict=True):
        direction = 0
        if "N" in (a, b):
            missing += 1
        elif a == b == "-":
            gap_gap += 1
        elif a == "-" or b == "-":
            gap_base += 1
            direction = 1 if a == "-" else 2
            tracts += direction != previous
        else:
            sites += 1
            snps += a != b
        previous = direction
    return snps, sites, gap_base, gap_gap, missing, tracts


@dataclass
class LocusAlignment:
    """Persisted, validated reconstruction with pairwise evidence arithmetic."""

    samples: dict
    partitions: list
    provenance: dict
    aligned_alleles: dict
    sample_ids: tuple = ()

    def pairs(self, min_loci=599):
        if not self.samples:
            raise ValueError("Pair iteration requires materialized samples")
        if not 1 <= min_loci <= len(self.partitions):
            raise ValueError("Minimum compared loci lies outside the canonical scheme")
        vectors = {ident: np.frombuffer(seq.encode(), dtype=np.uint8)
                   for ident, seq in self.samples.items()}
        for left, right in combinations(self.samples, 2):
            a, b = vectors[left], vectors[right]
            missing = (a == ord("N")) | (b == ord("N"))
            agap, bgap = a == ord("-"), b == ord("-")
            joint = ~(missing | agap | bgap)
            gap_base = ~missing & (agap ^ bgap)
            gap_gap = ~missing & agap & bgap
            directions = np.zeros(len(a), dtype=np.uint8)
            directions[gap_base & agap] = 1
            directions[gap_base & bgap] = 2
            compared = tracts = 0
            oracle = np.zeros(6, dtype=np.int64)
            for partition in self.partitions:
                if not partition["width"]:
                    continue
                start, end = partition["start"] - 1, partition["end"]
                compared += bool(np.any(joint[start:end]))
                part = directions[start:end]
                tracts += int(np.count_nonzero((part != 0) & (part != np.r_[0, part[:-1]])))
                oracle += _scalar_stats(self.samples[left][start:end], self.samples[right][start:end])
            counts = (int(np.count_nonzero((a != b) & joint)), int(joint.sum()),
                      int(gap_base.sum()), int(gap_gap.sum()), int(missing.sum()), tracts)
            if tuple(oracle) != counts or sum(counts[1:5]) != len(a):
                raise ValueError("Independent per-locus and concatenated arithmetic disagree")
            yield {"sample_a": left, "sample_b": right, "locus_snps": counts[0],
                   "compared_loci": compared, "compared_sites": counts[1],
                   "eligible": compared >= min_loci, "gap_vs_base_columns": counts[2],
                   "both_gap_columns": counts[3], "missing_columns": counts[4],
                   "indel_tracts": counts[5], "alignment_sites": len(a),
                   "locus_mismatch_fraction": counts[0] / counts[1] if counts[1] else None,
                   "independent_arithmetic_verified": True}


def _sample_profiles_hash(records, loci):
    """Hash the canonical JSON list without materializing all profile dictionaries."""
    digest = hashlib.sha256(b"[")
    for index, row in enumerate(records):
        if index:
            digest.update(b",")
        profile = _profile(row) or {}
        value = {"sample_id": row["sample_id"],
                 "calls": {locus: _allele(profile.get(locus)) for locus in loci}}
        digest.update(json.dumps(value, sort_keys=True, separators=(",", ":")).encode())
    digest.update(b"]")
    return digest.hexdigest()


def align_catalogue(catalogue, records, output, mafft, *, materialize_samples=True,
                    export_concatenation=True, progress=None):
    """Run one deterministic MSA per variable locus and reconstruct all samples.

    ``records`` must be canonical dataset-adapter rows in the desired sample
    order. Locus order is their explicit shared ``cgmlst_loci`` order. A locus
    without observed DNA retains a zero-width absent partition (null bounds).
    """
    sequences, scope, catalogue_provenance = read_catalogue(catalogue, records)
    return align_sequences(sequences, records, output, mafft, scope=scope,
                           provenance=catalogue_provenance,
                           materialize_samples=materialize_samples,
                           export_concatenation=export_concatenation, progress=progress)


def align_sequences(sequences, records, output, mafft, *, scope, provenance,
                    materialize_samples=True, export_concatenation=True, progress=None):
    """Align vetted ``(locus, allele) -> DNA`` with supplied source provenance.

    The caller is responsible for establishing sequence-source compatibility.
    This lower-level API records the supplied scope and provenance verbatim;
    catalogue/CLI consumers must continue through ``align_catalogue`` validation.
    ``materialize_samples=False`` keeps only allele alignments in memory and
    streams the optional concatenated FASTA one sample at a time.
    """
    catalogue_provenance = provenance
    if not isinstance(scope, dict) or not isinstance(provenance, dict) or not provenance:
        raise ValueError("Vetted alignments require explicit source scope and provenance")
    if any(not isinstance(dna, str) or not dna or set(dna) - set("ACGT")
           for dna in sequences.values()):
        raise ValueError("Vetted sequences must contain nonempty uppercase ACGT DNA")
    if not records or len({r["sample_id"] for r in records}) != len(records):
        raise ValueError("Require distinct prepared samples")
    loci = list(records[0]["cgmlst_loci"])
    if not loci or len(set(loci)) != len(loci) or any(list(r["cgmlst_loci"]) != loci for r in records):
        raise ValueError("Samples must share an explicit ordered canonical locus universe")
    for row in records:
        if any(c.isspace() for c in row["sample_id"]) or not row["sample_id"]:
            raise ValueError("Sample identifiers must be nonempty FASTA identifiers")
    executable = Path(mafft).expanduser().resolve()
    try:
        version_run = subprocess.run([str(executable), "--version"], capture_output=True, text=True,
                                     check=True, timeout=30)
    except (OSError, subprocess.SubprocessError) as error:
        raise ValueError(f"MAFFT version check failed: {error}") from error
    version = (version_run.stdout + version_run.stderr).strip()
    if not version:
        raise ValueError("MAFFT returned no version provenance")
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    # A failed rerun must never leave a complete marker beside partial new files.
    (output / "alignment.json").unlink(missing_ok=True)
    chunks = {row["sample_id"]: [] for row in records} if materialize_samples else None
    partitions, missing_calls, commands, aligned_alleles = [], [], [], {}
    offset = 0
    for number, locus in enumerate(loci, 1):
        directory = output / "loci" / f"{number:04d}"
        directory.mkdir(parents=True, exist_ok=True)
        observed = {}
        for row in records:
            allele = _allele((_profile(row) or {}).get(locus))
            if allele is not None and (locus, allele) in sequences:
                observed[allele] = sequences[locus, allele]
        unique = sorted(set(observed.values()))
        originals = {f"dna_{i:04d}": dna for i, dna in enumerate(unique, 1)}
        sequence_ids = {dna: ident for ident, dna in originals.items()}
        mapping = {allele: sequence_ids[dna] for allele, dna in sorted(observed.items())}
        save_json(directory / "allele-map.json", mapping)
        input_path, aligned_path = directory / "input.fasta", directory / "aligned.fasta"
        input_path.write_text(fasta_text(originals))
        if len(unique) > 1:
            args = [str(executable), *MAFFT_OPTIONS, str(input_path)]
            try:
                run = subprocess.run(args, capture_output=True, text=True, check=True, timeout=180)
            except (OSError, subprocess.SubprocessError) as error:
                raise ValueError(f"MAFFT failed for canonical locus {locus}: {error}") from error
            aligned = read_fasta(run.stdout)
            width = validate_msa(aligned, originals)
            commands.append({"locus": locus, "args": args, "stderr": run.stderr})
            method = "mafft"
        elif unique:
            aligned, width, method = originals, len(unique[0]), "identical-sequence-exact"
            validate_msa(aligned, originals)
        else:
            aligned, width, method = {}, 0, "absent-no-observed-dna"
        aligned_path.write_text(fasta_text(aligned))
        aligned_alleles.update({(locus, allele): aligned[identifier] for allele, identifier in mapping.items()})
        partition = {"locus": locus, "ordinal": number, "start": offset + 1 if width else None,
                     "end": offset + width if width else None, "width": width,
                     "status": "represented" if width else "absent", "method": method,
                     "observed_alleles": len(observed), "unique_sequences": len(unique),
                     "input_path": str(input_path.relative_to(output)), "input_sha256": file_hash(input_path),
                     "aligned_path": str(aligned_path.relative_to(output)), "aligned_sha256": file_hash(aligned_path),
                     "allele_map_sha256": file_hash(directory / "allele-map.json")}
        partitions.append(partition)
        for row in records:
            ident = row["sample_id"]
            allele = _allele((_profile(row) or {}).get(locus))
            if allele in mapping:
                sequence = aligned[mapping[allele]]
            else:
                sequence = "N" * width
                missing_calls.append({"sample_id": ident, "locus": locus, "allele": allele,
                                      "reason": "missing_call" if allele is None else "missing_sequence",
                                      "alignment_width": width})
            if chunks is not None:
                chunks[ident].append(sequence)
        offset += width
        if progress and (number % 25 == 0 or number == len(loci)):
            progress("align_loci", number, len(loci))
    samples = {ident: "".join(parts) for ident, parts in chunks.items()} if chunks else {}
    if not offset or any(len(sequence) != offset for sequence in samples.values()):
        raise ValueError("Reconstructed sample alignment is empty or inconsistent")
    concat = output / "concatenated.fasta"
    if export_concatenation:
        with concat.open("w") as stream:
            for row in records:
                ident = row["sample_id"]
                stream.write(f">{ident}\n")
                written = 0
                for partition in partitions:
                    locus, width = partition["locus"], partition["width"]
                    allele = _allele((_profile(row) or {}).get(locus))
                    sequence = aligned_alleles.get((locus, allele), "N" * width)
                    stream.write(sequence)
                    written += len(sequence)
                stream.write("\n")
                if written != offset:
                    raise ValueError("Persisted sample alignment has inconsistent width")
        if materialize_samples and read_fasta(concat.read_text()) != samples:
            raise ValueError("Persisted concatenation differs from reconstructed samples")
    else:
        concat.unlink(missing_ok=True)
    save_json(output / "partitions.json", partitions)
    save_json(output / "missing-calls.json", missing_calls)
    provenance = {"schema": "chronoclade.cgmlst.locus-alignment", "schema_version": 1,
                  "status": "complete", "scope": scope, "catalogue": catalogue_provenance,
                  "implementation_sha256": file_hash(__file__),
                  "sample_profiles_sha256": _sample_profiles_hash(records, loci),
                  "mafft": {"path": str(executable), "sha256": file_hash(executable), "version": version,
                            "options": list(MAFFT_OPTIONS), "commands": commands},
                  "canonical_loci": len(loci), "aligned_variable_loci": len(commands),
                  "alignment_sites": offset, "sample_ids": [row["sample_id"] for row in records],
                  "materialized_samples": materialize_samples,
                  "concatenation_exported": export_concatenation,
                  "concatenation_sha256": file_hash(concat) if export_concatenation else None,
                  "partitions_sha256": file_hash(output / "partitions.json"),
                  "missing_calls_sha256": file_hash(output / "missing-calls.json"),
                  "missing_policy": "N across full locus alignment; globally absent loci have zero width and null coordinate bounds",
                  "pair_policy": "Substitutions only at jointly ACGT columns; gaps and missingness reported separately; indel tracts reset at locus boundaries and when gap direction changes"}
    save_json(output / "alignment.json", provenance)
    return LocusAlignment(samples, partitions, provenance, aligned_alleles,
                          tuple(row["sample_id"] for row in records))


def load_alignment(output, records, *, materialize_samples=False):
    """Resume a completed alignment after checking persisted evidence hashes.

    Every per-locus input, alignment and allele map is verified and the supplied
    sample profiles must reproduce the recorded order and fingerprint. Large
    concatenations are hashed incrementally and loaded only when requested.
    Recorded implementation provenance is preserved across software updates.
    """
    output = Path(output).resolve()
    provenance = json.loads((output / "alignment.json").read_text())
    if (provenance.get("schema") != "chronoclade.cgmlst.locus-alignment"
            or provenance.get("schema_version") != 1 or provenance.get("status") != "complete"):
        raise ValueError("Require a completed supported locus alignment")
    ids = tuple(row["sample_id"] for row in records)
    if list(ids) != provenance.get("sample_ids") or len(set(ids)) != len(ids):
        raise ValueError("Persisted alignment sample order differs from supplied profiles")

    def verified_path(relative, expected):
        path = (output / relative).resolve()
        if not path.is_relative_to(output) or file_hash(path) != expected:
            raise ValueError(f"Persisted alignment artifact hash mismatch: {relative}")
        return path

    partitions_path = verified_path("partitions.json", provenance["partitions_sha256"])
    verified_path("missing-calls.json", provenance["missing_calls_sha256"])
    partitions = json.loads(partitions_path.read_text())
    loci = [partition["locus"] for partition in partitions]
    if (not loci or len(set(loci)) != len(loci)
            or any(list(row["cgmlst_loci"]) != loci for row in records)
            or len(loci) != provenance["canonical_loci"]
            or _sample_profiles_hash(records, loci) != provenance["sample_profiles_sha256"]):
        raise ValueError("Persisted alignment canonical loci or sample profiles differ")
    aligned_alleles, offset = {}, 0
    for partition in partitions:
        input_path = verified_path(partition["input_path"], partition["input_sha256"])
        aligned_path = verified_path(partition["aligned_path"], partition["aligned_sha256"])
        map_path = verified_path(str(Path(partition["input_path"]).parent / "allele-map.json"),
                                 partition["allele_map_sha256"])
        original, aligned = read_fasta(input_path.read_text()), read_fasta(aligned_path.read_text())
        mapping = json.loads(map_path.read_text())
        width = validate_msa(aligned, original) if original else 0
        if (width != partition["width"] or (not width and (aligned or mapping))
                or partition["start"] != (offset + 1 if width else None)
                or partition["end"] != (offset + width if width else None)
                or any(identifier not in aligned for identifier in mapping.values())):
            raise ValueError("Persisted per-locus alignment and partition evidence disagree")
        aligned_alleles.update({(partition["locus"], allele): aligned[identifier]
                               for allele, identifier in mapping.items()})
        offset += width
    if offset != provenance["alignment_sites"]:
        raise ValueError("Persisted alignment width differs from provenance")
    samples = {}
    if provenance.get("concatenation_exported", True):
        concat = verified_path("concatenated.fasta", provenance["concatenation_sha256"])
        if materialize_samples:
            samples = read_fasta(concat.read_text())
            if tuple(samples) != ids or any(len(sequence) != offset for sequence in samples.values()):
                raise ValueError("Persisted concatenation sample order or width differs")
    elif materialize_samples:
        for row in records:
            profile = _profile(row) or {}
            samples[row["sample_id"]] = "".join(
                aligned_alleles.get((part["locus"], _allele(profile.get(part["locus"]))),
                                    "N" * part["width"]) for part in partitions)
    return LocusAlignment(samples, partitions, provenance, aligned_alleles, ids)
