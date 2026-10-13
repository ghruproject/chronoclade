"""Scope-checked nucleotide comparisons of frozen cgMLST allele sequences.

Only equal-length alleles are compared sitewise. Missing DNA and length-changing
alleles are unavailable evidence, never inferred from allele identifiers.
"""

from __future__ import annotations

from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path
import re

import numpy as np

from chronoclade.artifacts import write_json
from chronoclade.context_refinement import _allele, _loci, _profile
from chronoclade.profile_analysis import _pair

CATALOGUE_SCHEMA = "chronoclade.allele-sequence-catalogue"
_SCOPE_FIELDS = ("scheme_id", "scheme_version", "database_version", "database_sha256")


def read_catalogue(path, records):
    """Validate bytes and exact species/scheme/database scope before analysis."""
    path = Path(path).resolve()
    try:
        raw = path.read_bytes()
        manifest = json.loads(raw)
        if (manifest.get("schema") != CATALOGUE_SCHEMA
                or type(manifest.get("schema_version")) is not int
                or manifest["schema_version"] != 1 or manifest.get("status") != "complete"):
            raise ValueError("Expected a complete chronoclade.allele-sequence-catalogue v1")
        scope = manifest["scope"]
        if not isinstance(scope, dict) or set(scope) != {"species", *_SCOPE_FIELDS}:
            raise ValueError("Sequence catalogue requires explicit species and database scope")
        for key in ("species", "scheme_id", "scheme_version"):
            if not isinstance(scope[key], str) or not scope[key].strip():
                raise ValueError(f"Sequence catalogue requires nonempty {key}")
        if scope["database_version"] is not None and (
                not isinstance(scope["database_version"], str) or not scope["database_version"].strip()):
            raise ValueError("Catalogue database_version must be text or null")
        if not isinstance(scope["database_sha256"], str) or not re.fullmatch(
                r"[0-9a-f]{64}", scope["database_sha256"]):
            raise ValueError("Sequence catalogue requires a known database SHA256 fingerprint")
        loci = set()
        for row in records:
            if row.get("species") != scope["species"]:
                raise ValueError("Sequence catalogue species is incompatible with prepared samples")
            if _profile(row) is None:
                continue
            for field in _SCOPE_FIELDS:
                if row.get("cgmlst_" + field.removesuffix("_id")) != scope[field]:
                    raise ValueError(f"Sequence catalogue {field} is incompatible with prepared data")
            loci.update(row.get("cgmlst_loci") or [])
        descriptor = manifest["artifact"]
        relative = Path(descriptor["path"])
        if relative.is_absolute() or ".." in relative.parts or not relative.parts:
            raise ValueError("Catalogue DNA CSV path must be a local relative path")
        dna_path = (path.parent / relative).resolve()
        if not dna_path.is_relative_to(path.parent):
            raise ValueError("Catalogue DNA CSV must remain within its manifest directory")
        csv_bytes = dna_path.read_bytes()
        digest = hashlib.sha256(csv_bytes).hexdigest()
        if descriptor["sha256"] != digest:
            raise ValueError("Catalogue DNA CSV SHA256 mismatch")
        reader = csv.DictReader(io.StringIO(csv_bytes.decode("utf-8"), newline=""))
        if reader.fieldnames != ["locus", "allele", "dna"]:
            raise ValueError("Catalogue CSV requires exactly locus,allele,dna columns in order")
        sequences = {}
        for number, row in enumerate(reader, 2):
            if None in row or any(not isinstance(row.get(key), str) or not row[key]
                                  or row[key] != row[key].strip() for key in reader.fieldnames):
                raise ValueError(f"Catalogue CSV row {number} has missing or malformed values")
            locus, allele, dna = row["locus"], row["allele"], row["dna"].upper()
            if locus not in loci:
                raise ValueError(f"Catalogue locus {locus!r} is outside the prepared scheme")
            if _allele(allele) is None:
                raise ValueError(f"Catalogue row {number} uses a missing marker as an allele")
            if set(dna) - set("ACGT"):
                raise ValueError(f"Catalogue row {number} has ambiguous or unsupported DNA bases")
            key = (locus, allele)
            if key in sequences:
                raise ValueError(f"Catalogue repeats locus/allele {key!r}")
            sequences[key] = dna
        return sequences, scope, {
            "manifest_path": str(path), "manifest_sha256": hashlib.sha256(raw).hexdigest(),
            "dna_csv_path": str(dna_path), "dna_csv_sha256": digest,
            "catalogue_records": len(sequences),
        }
    except (OSError, UnicodeError, KeyError, TypeError, AttributeError, csv.Error) as error:
        raise ValueError(f"Invalid frozen allele-sequence catalogue: {error}") from error


class DNAComparator:
    """A pair/locus-resampling adapter for the existing descriptive profile engine."""

    metadata = {
        "distance_method": "direct-dna",
        "distance_definition": "fraction of nucleotide substitutions across jointly called, equal-length cgMLST allele sequences",
        "distance_units": "cgMLST nucleotide mismatch fraction",
        "raw_difference_field": "snp_differences",
        "nearest_ranking": "raw nucleotide substitutions among sufficiently jointly callable, equal-length cgMLST alleles; all raw-distance ties retained",
    }
    limitations = [
        "Direct DNA compares equal-length cgMLST alleles sitewise, assuming homologous positions; no alignment is inferred.",
        "Missing DNA and unequal-length allele pairs are excluded; every comparison records its loci and nucleotide denominator.",
        "These are substitutions within supplied cgMLST genes, not genome-wide or recombination-corrected SNP distances.",
    ]

    def __init__(self, path, records):
        self.sequences, self.scope, self.provenance = read_catalogue(path, records)
        self._locus_cache, self._pair_cache = {}, {}
        self.audit = []

    def locus_pair(self, locus, left, right):
        key = (locus, left, right)
        if key not in self._locus_cache:
            if left is None or right is None:
                result = (0, 0, "missing_call")
            else:
                a, b = self.sequences.get((locus, left)), self.sequences.get((locus, right))
                if a is None or b is None:
                    result = (0, 0, "missing_sequence")
                elif len(a) != len(b):
                    result = (0, 0, "unequal_length")
                else:
                    result = (sum(x != y for x, y in zip(a, b, strict=True)), len(a), None)
            self._locus_cache[key] = result
        return self._locus_cache[key]

    def pair(self, left, right, overlap):
        key = (left["sample_id"], right["sample_id"], overlap)
        if key in self._pair_cache:
            return self._pair_cache[key]
        evidence, reason = _pair(left, right, overlap)
        if evidence is None:
            return None, reason
        a, b = _profile(left), _profile(right)
        loci = sorted(_loci(left, a))
        failures = {reason: [] for reason in ("missing_call", "missing_sequence", "unequal_length")}
        differences = sites = compared = alleles = 0
        for locus in loci:
            x, y = _allele(a.get(locus)), _allele(b.get(locus))
            difference, count, reason = self.locus_pair(locus, x, y)
            if reason:
                failures[reason].append(locus)
            else:
                differences += difference
                sites += count
                compared += 1
                alleles += x != y
        eligible = bool(compared and compared / len(loci) >= overlap)
        result = {
            "allele_differences": alleles, "snp_differences": differences,
            "shared_called_sites": sites, "shared_called_loci": compared,
            "scheme_loci": len(loci), "call_overlap": compared / len(loci),
            "distance": differences / sites if sites else None,
            "missing_call_loci": failures["missing_call"],
            "missing_sequence_loci": failures["missing_sequence"],
            "unequal_length_loci": failures["unequal_length"],
        }
        if any(failures.values()):
            self.audit.append({"sample_a": key[0], "sample_b": key[1],
                               "status": "comparable" if eligible else "unavailable", **result})
        answer = (result, None) if eligible else (None, "insufficient_sequence_overlap")
        self._pair_cache[key] = answer
        return answer

    def resample(self, rows, loci, selected, overlap):
        """Resample whole genes with O(samples squared) working memory."""
        profiles = [_profile(row) for row in rows]
        n = len(rows)
        differences = np.zeros((n, n))
        sites = np.zeros_like(differences)
        calls = np.zeros_like(differences)
        for index, multiplicity in Counter(selected).items():
            locus = loci[index]
            categories = list(dict.fromkeys(_allele(profile.get(locus)) for profile in profiles))
            positions = {allele: i for i, allele in enumerate(categories)}
            encoded = np.array([positions[_allele(profile.get(locus))] for profile in profiles])
            shape = (len(categories), len(categories))
            raw, bases, present = (np.zeros(shape) for _ in range(3))
            for i, left in enumerate(categories):
                for j in range(i + 1):
                    diff, count, reason = self.locus_pair(locus, left, categories[j])
                    for array, value in ((raw, diff), (bases, count), (present, int(reason is None))):
                        array[i, j] = array[j, i] = value
            positions = np.ix_(encoded, encoded)
            differences += multiplicity * raw[positions]
            sites += multiplicity * bases[positions]
            calls += multiplicity * present[positions]
        if not np.all(calls / len(selected) >= overlap) or np.any(sites <= 0):
            return None
        return differences / sites

    def save_audit(self, directory):
        path = Path(directory) / "sequence_comparison_audit.json"
        write_json(path, {"schema": "chronoclade.cgmlst.sequence-comparison-audit", "schema_version": 1,
                          "scope": self.scope, "provenance": self.provenance,
                          "pair_exclusions": self.audit,
                          "exclusion_counts": dict(Counter(
                              reason for row in self.audit for reason in ("missing_call_loci", "missing_sequence_loci", "unequal_length_loci")
                              if row[reason]))})
        return path
