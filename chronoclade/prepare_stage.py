"""Independent preparation of frozen or live query-only typed datasets."""

from __future__ import annotations

import json
import re
import shutil
import tempfile
import shlex
import uuid
from dataclasses import dataclass, replace
from html import escape
from pathlib import Path

from chronoclade import __version__
from chronoclade.artifacts import file_sha256, write_json
from chronoclade.datasets import (
    SCHEMA_NAME,
    DatasetError,
    LocusCatalogue,
    from_profile_records,
    load_dataset,
    write_dataset,
)
from chronoclade.datasets.model import EVIDENCE_TABLES, check_json
from chronoclade.report_components.styles import report_styles
from chronoclade.typing_scopes import typing_scope


@dataclass(frozen=True)
class PrepareResult:
    dataset_manifest: Path
    report_path: Path
    audit_path: Path
    stage_manifest: Path


def _json(path: Path):
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
        check_json(value, location=path.name)
        return value
    except (OSError, ValueError) as error:
        raise DatasetError(f"Cannot read preparation JSON {path}: {error}") from error


def _merge_frozen_metadata(dataset, path):
    """Apply exact-ID metadata, including supplied classical ST annotations."""
    from chronoclade.datasets.adapters import samples_from_records
    from chronoclade.profile_inputs import _csv

    allowed = {
        "country", "region", "city", "nuts2", "host", "isolation_source", "collection_date",
        "hospital", "hospital_code", "publication_isolate_id", "specimen_type",
        "clinical_specimen", "study_accession", "patient_id", "location",
        "mlst_st", "mlst_scheme",
    }
    identity = {"accession", "source_genome_id", "sample_id"}
    crosswalk = {r["sample_id"]: r for r in dataset.crosswalk}
    samples = [dict(r) for r in dataset.samples]
    provenance = [dict(r) for r in dataset.provenance]
    conflicts = list(dataset.conflicts)
    for override in _csv(path):
        unknown = set(override) - allowed - identity
        if unknown:
            raise DatasetError(f"Unsupported metadata override fields: {sorted(unknown)}")
        identifiers = {override[k] for k in identity if override.get(k)}
        identifier = ", ".join(sorted(identifiers)) or None
        matches = []
        for row in samples:
            evidence = crosswalk.get(row["sample_id"], {})
            aliases = {row["sample_id"], row.get("source_genome_id"), evidence.get("source_genome_id")}
            for source in (row, evidence):
                for key in ("aliases", "run_accessions", "biosample_accessions", "assembly_accessions"):
                    aliases.update(source.get(key) or [])
                for key in ("accession", "biosample", "sample_accession", "run_accession", "assembly_accession"):
                    if source.get(key):
                        aliases.add(source[key])
            if identifiers and identifiers <= aliases:
                matches.append(row)
        if len(matches) != 1:
            raise DatasetError(f"Metadata override {identifier!r} does not match exactly one sample")
        row = matches[0]
        for field in sorted(allowed & set(override)):
            value = override[field]
            if not value:
                continue
            if field == 'mlst_st' and not re.fullmatch(r'[1-9][0-9]*', value):
                raise DatasetError('Metadata mlst_st must be a positive numeric ST')
            if row.get(field) not in (None, "", value):
                conflicts.append({"sample_id": row["sample_id"], "field": field,
                                  "source": "prepared_metadata", "value": row[field]})
            if field in {'mlst_st', 'mlst_scheme'}:
                # User annotations must not retain alleles/provenance from a
                # different earlier native call as evidence for the new ST.
                for stale in ('mlst_alleles', 'mlst_database_sha256', 'mlst_tool_version',
                              'mlst_assembly_sha256', 'mlst_reason'):
                    row.pop(stale, None)
                row.update(mlst_source='user_metadata', mlst_status='supplied')
            row[field] = value
            for previous in provenance:
                if previous.get("sample_id") == row["sample_id"] and previous.get("field") == field and previous.get("selected"):
                    previous.update(selected=False, superseded_by="user_metadata_override")
            provenance.append({"sample_id": row["sample_id"], "field": field,
                               "source": "user", "value": value, "selected": True,
                               "metadata_sha256": file_sha256(path)})
        if override.get("collection_date"):
            for field in ("date_start", "date_end", "date_precision", "date_raw"):
                row.pop(field, None)
            normalized = samples_from_records([row])[0]
            for field in ("date_start", "date_end", "date_precision", "date_raw"):
                row[field] = normalized[field]
    return replace(dataset, samples=tuple(samples), provenance=tuple(provenance),
                   conflicts=tuple(conflicts), dataset_id=None)


def _catalogues(path: Path) -> tuple[LocusCatalogue, ...]:
    raw = _json(path)
    if isinstance(raw, dict):
        if set(raw) != {"catalogues"}:
            raise DatasetError("Catalogue JSON object requires only a catalogues list")
        raw = raw["catalogues"]
    if not isinstance(raw, list) or not raw:
        raise DatasetError("Catalogue JSON requires a nonempty catalogue list")
    result = []
    for index, value in enumerate(raw):
        if not isinstance(value, dict) or not isinstance(value.get("loci"), list):
            raise DatasetError(f"catalogues[{index}] requires an explicit loci list")
        if "complete" not in value:
            raise DatasetError(f"catalogues[{index}] requires explicit complete true/false")
        try:
            catalogue = LocusCatalogue(**dict(value, loci=tuple(value["loci"])))
            catalogue.validate()
        except (TypeError, ValueError) as error:
            raise DatasetError(f"catalogues[{index}]: {error}") from error
        result.append(catalogue)
    return tuple(result)


def _frozen_dataset(raw, catalogues, species):
    evidence = {}
    if isinstance(raw, dict):
        allowed = {"records", "parameters", *EVIDENCE_TABLES}
        if set(raw) - allowed or "records" not in raw:
            raise DatasetError(
                "Frozen input requires records and optional evidence tables/parameters; "
                "legacy queries/context resolver outputs need explicit adaptation"
            )
        evidence = {key: value for key, value in raw.items() if key != "records"}
        raw = raw["records"]
    if not isinstance(raw, list) or not raw or any(not isinstance(row, dict) for row in raw):
        raise DatasetError("Frozen input requires a nonempty list of canonical record objects")
    for name in EVIDENCE_TABLES:
        if name in evidence and not isinstance(evidence[name], list):
            raise DatasetError(f"Frozen {name} must be an explicit list")
    if "parameters" in evidence and not isinstance(evidence["parameters"], dict):
        raise DatasetError("Frozen parameters must be an object")
    records = []
    for index, row in enumerate(raw):
        if not isinstance(row.get("sample_id"), str):
            raise DatasetError(f"records[{index}] requires a stable sample_id")
        role = row.get("role")
        origin = row.get("origin")
        if role is None and origin is None:
            role = "input"
        record = dict(row, role=role)
        if species:
            existing = row.get("species")
            if existing and (
                not isinstance(existing, str) or existing.casefold() != species.casefold()
            ):
                raise DatasetError(f"Sample {row['sample_id']!r} species differs from --species")
            record["species"] = species
        records.append(record)
    try:
        return from_profile_records(records, catalogues, **evidence)
    except (KeyError, TypeError, AttributeError) as error:
        raise DatasetError(f"Malformed frozen profile evidence: {error}") from error


def _portable_assemblies(dataset, source_root, staging):
    """Copy explicitly linked local assemblies; leave external identities unresolved."""
    rows, artifacts = [], {}
    for row in dataset.samples:
        row = dict(row)
        reference = row.get("assembly_reference")
        if reference:
            if "://" in reference:
                rows.append(row)
                continue
            path = Path(reference).expanduser()
            path = path if path.is_absolute() else source_root / path
            if path.is_file():
                digest = file_sha256(path)
                relative = f"assemblies/{digest}{path.suffix}"
                destination = staging / relative
                destination.parent.mkdir(exist_ok=True)
                if not destination.exists():
                    shutil.copyfile(path, destination)
                row["assembly_reference"] = relative
                artifacts[relative] = {
                    "path": relative,
                    "sha256": digest,
                    "bytes": destination.stat().st_size,
                }
            elif Path(reference).is_absolute() or "/" in reference or Path(reference).suffix:
                raise DatasetError(
                    f"Sample {row['sample_id']!r} assembly file is missing: {reference}"
                )
        rows.append(row)
    return replace(dataset, samples=tuple(rows)), list(artifacts.values())


def _validate_previous_stage(root):
    """Validate extra assembly/source/report artifacts when importing our output."""
    root = root.resolve()
    path = root / "prepare.json"
    if not path.exists():
        return
    stage = _json(path)
    if (
        not isinstance(stage, dict)
        or stage.get("schema") != "chronoclade.prepare-result"
        or type(stage.get("schema_version")) is not int
        or stage["schema_version"] != 1
        or stage.get("status") != "complete"
        or not isinstance(stage.get("artifacts"), list)
    ):
        raise DatasetError("Existing prepare result requires a complete supported manifest")
    seen = set()
    for item in stage["artifacts"]:
        if not isinstance(item, dict):
            raise DatasetError("Existing prepare artifact requires a descriptor")
        raw = item.get("path")
        if (
            not isinstance(raw, str)
            or not raw
            or "\\" in raw
            or Path(raw).is_absolute()
            or any(part in {".", "..", ""} for part in raw.split("/"))
        ):
            raise DatasetError("Existing prepare artifact requires a safe relative path")
        artifact = (root / raw).resolve()
        if not artifact.is_relative_to(root) or raw in seen or not artifact.is_file():
            raise DatasetError(
                "Existing prepare artifact is missing, duplicated or outside the bundle"
            )
        seen.add(raw)
        if (
            type(item.get("bytes")) is not int
            or artifact.stat().st_size != item["bytes"]
            or file_sha256(artifact) != item.get("sha256")
        ):
            raise DatasetError(f"Existing prepare artifact checksum/size mismatch: {raw}")
    if not {"dataset.json", "readiness.json", "report.html"}.issubset(seen):
        raise DatasetError("Existing prepare result is missing its required artifacts")


def _expected_lineage(species):
    species = (species or "").casefold()
    if "klebsiella" in species:
        return "cglin"
    if species in {"escherichia coli", "e. coli", "ecoli"}:
        return "hiercc"
    return None


def _lineage_ready(evidence, kind):
    """Retain the provider's complete cgLIN status without accepting provisional codes."""
    status = evidence.get("resolution_status")
    if status not in ({"resolved", "complete"} if kind == "cglin" else {"resolved"}):
        return False
    frozen = evidence.get("evidence") or {}
    if not isinstance(frozen, dict):
        raise DatasetError(f"Sample {evidence['sample_id']!r} lineage evidence must be an object")
    provisional = frozen.get(kind + "_provisional")
    if provisional is True or str(provisional).casefold() in {"true", "1", "yes"}:
        return False
    scope = typing_scope(
        dict(
            frozen,
            **{
                kind + "_scheme": evidence.get("scheme_id"),
                kind + "_scheme_version": evidence.get("scheme_version"),
            },
        ),
        kind,
    )
    return bool(evidence.get("kind") == kind and evidence.get("assignment") and scope)


def _readiness(dataset, manifest):
    profiles, lineages = {}, {}
    for matrix in dataset.profiles:
        for ident, codes in zip(matrix.sample_ids, matrix.codes, strict=True):
            profiles.setdefault(ident, []).append(
                {
                    "scheme_id": matrix.catalogue.scheme_id,
                    "scheme_version": matrix.catalogue.scheme_version,
                    "database_version": matrix.catalogue.database_version,
                    "database_sha256": matrix.catalogue.database_sha256,
                    "catalogue_complete": matrix.catalogue.complete,
                    "loci": len(matrix.catalogue.loci),
                    "called_loci": int((codes != 0).sum()),
                    "missing_loci": int((codes == 0).sum()),
                }
            )
    for row in dataset.lineages:
        lineages.setdefault(row["sample_id"], []).append(row)
    samples = []
    for row in dataset.samples:
        ident = row["sample_id"]
        profile_rows = profiles.get(ident, [])
        kind = _expected_lineage(row.get("species"))
        lineage_rows = lineages.get(ident, [])
        lineage_available = (
            any(_lineage_ready(evidence, kind) for evidence in lineage_rows) if kind else False
        )
        mlst_st = str(row.get("mlst_st") or "").strip()
        mlst_st_ready = bool(re.fullmatch(r"[1-9][0-9]*", mlst_st))
        reasons = []
        if not any(evidence["called_loci"] for evidence in profile_rows):
            reasons.append(
                "cgMLST profile unavailable after configured profile calling"
                if str(
                    (dataset.parameters.get("prepare_operations") or {}).get("typing", "")
                ).startswith("run")
                else "cgMLST profile unavailable; no profile calling was performed"
            )
        if any(not evidence["catalogue_complete"] for evidence in profile_rows):
            reasons.append("Locus catalogue is explicitly incomplete")
        if kind and not lineage_available:
            reasons.append(f"Resolved {kind} assignment with a known scope is unavailable")
        if not kind:
            reasons.append("Lineage readiness is not assessed for this species")
        if row['role'] == 'input' and not mlst_st_ready:
            reasons.append(
                "A verified seven-gene MLST ST is required for public context discovery; "
                "prepare calls classical MLST for missing STs when an input assembly and mapped caller/database are available. "
                f"Current status: {row.get('mlst_reason') or 'not available'}. "
                "Supply an assembly and install the caller with `pixi install --frozen -e mlst`, "
                "or add numeric `mlst_st` from verified metadata before `cgmlst --fetch-context`. "
                "Do not infer it from cgLIN/HierCC or nearest-neighbour labels."
            )
        reference = row.get("assembly_reference")
        samples.append(
            {
                "sample_id": ident,
                "label": row["label"],
                "role": row["role"],
                "profiles": profile_rows,
                "required_lineage_kind": kind,
                "lineage_ready": bool(lineage_available),
                "mlst_st_ready": mlst_st_ready,
                "mlst_st": row.get('mlst_st'),
                "mlst_scheme": row.get('mlst_scheme'),
                "mlst_status": row.get('mlst_status') or ('supplied' if mlst_st_ready else 'unavailable'),
                "profile_ready": bool(any(evidence["called_loci"] for evidence in profile_rows)),
                "readiness": "ready" if not reasons else "incomplete",
                "reasons": reasons,
                "assembly_status": "bundled"
                if reference and (manifest.parent / reference).is_file()
                else "external_reference"
                if reference
                else "retrieval_not_assessed"
                if row.get("source_genome_id") or row.get("numeric_source_id")
                else "unavailable",
                "missing_metadata": [
                    field
                    for field in ("country", "region", "city", "nuts2", "host", "isolation_source")
                    if not row.get(field)
                ]
                + (
                    ["collection_date"]
                    if row.get("date_precision") in {None, "missing", "invalid"}
                    else []
                ),
            }
        )
    return {
        "schema": "chronoclade.prepare-readiness",
        "schema_version": 1,
        "status": "complete",
        "dataset_manifest": "dataset.json",
        "dataset_sha256": file_sha256(manifest),
        "dataset_id": json.loads(manifest.read_text())["dataset_id"],
        "sample_ids": list(dataset.sample_ids),
        "input_ids": [row["sample_id"] for row in dataset.samples if row["role"] == "input"],
        "context_ids": [row["sample_id"] for row in dataset.samples if row["role"] == "context"],
        "readiness": "ready"
        if all(row["readiness"] == "ready" for row in samples)
        else "incomplete",
        "samples": samples,
        "provider_context": dataset.parameters.get("profile_context"),
        "operations": dataset.parameters.get("prepare_operations")
        or {
            "typing": "not_run_frozen_import",
            "lineage_assignment": "not_run_frozen_import",
            "provider_retrieval": "not_run_offline",
            "context_discovery": "not_run",
            "typing_tool_installation": "not_checked_not_required_for_import",
            "reference_database": "not_checked_not_required_for_import",
        },
    }


def _operation_label(value):
    """Translate operation evidence for readers without changing saved statuses."""
    return {
        "run_exact_queries": "Exact query retrieval attempted",
        "not_run_local_assemblies": "Not requested; local assemblies supplied",
        "not_run_offline": "Not requested; frozen evidence imported",
        "run_missing_profile_or_lineage_inputs": "Native typing run where needed",
        "run_missing_profile_inputs": "Called missing input cgMLST profiles",
        "not_run_frozen_import": "Not run; frozen typing imported",
        "not_run": "Not run",
        "ready_validated_configuration": "Configured and validated",
        "unavailable_not_configured": "Not configured",
        "unavailable_configured_database_or_tool_not_ready": (
            "Unavailable; configured database or tool is not ready"
        ),
        "not_checked_not_required_for_import": (
            "Not checked; not required for imported evidence"
        ),
    }.get(value, value)


def _parallel_prepare_report(audit, parallel_prepare, evidence_root):
    """Explain worker results while leaving the collector's import status intact."""
    if not isinstance(parallel_prepare, dict) or not parallel_prepare.get('tasks'):
        return ''
    tasks = parallel_prepare['tasks']
    by_id = {row['sample_id']: row for row in audit['samples'] if row['role'] == 'input'}
    count = len(tasks)
    profiles = sum(by_id.get(task.get('sample_id'), {}).get('profile_ready', False) for task in tasks)
    sts = sum(by_id.get(task.get('sample_id'), {}).get('mlst_st_ready', False) for task in tasks)
    summary = (
        '<h3>Parallel preparation</h3>'
        f'<p>{count} checked per-assembly worker results were collected. '
        f'In the combined dataset, {profiles} / {count} inputs have cgMLST profiles '
        f'and {sts} / {count} have classical MLST STs. '
        'The collector imported validated worker results; it did not perform typing. '
        'The import statuses below describe collection, while worker evidence records '
        'the original calls or reused typing.</p>'
    )
    if evidence_root is None:
        return summary
    root = Path(evidence_root).resolve()
    entries = []
    for task in tasks:
        raw = task.get('retained_evidence')
        if not isinstance(raw, str) or not raw:
            continue
        directory = (root / raw).resolve()
        if not directory.is_relative_to(root):
            continue
        links = []
        for relative, label in (
            ('provider/native_typing/query_typing.json', 'Native cgMLST worker evidence'),
            ('classical_mlst/audit.json', 'Classical MLST worker evidence'),
        ):
            path = (directory / relative).resolve()
            if path.is_relative_to(root) and path.is_file():
                from urllib.parse import quote
                href = quote(path.relative_to(root).as_posix(), safe='/')
                links.append(f'<a href="{escape(href, quote=True)}">{label}</a>')
        if links:
            row = by_id.get(task.get('sample_id'), {})
            label = row.get('label') or task.get('sample_id') or 'Worker'
            entries.append(f'<li>{escape(str(label))}: ' + ' · '.join(links) + '</li>')
    if entries:
        summary += ('<details class="coverage-details"><summary>Worker typing evidence</summary>'
                    '<ul>' + ''.join(entries) + '</ul></details>')
    return summary


def _report(audit, *, parallel_prepare=None, evidence_root=None):
    def text(value):
        return escape(str(value), quote=True)

    ready = sum(row['profile_ready'] for row in audit['samples'])
    parallel_summary = _parallel_prepare_report(audit, parallel_prepare, evidence_root)
    typing_label = 'Collector cgMLST typing' if parallel_summary else 'cgMLST typing'
    eligibility = _next_stage_summary(audit)
    context_audit = audit.get("provider_context")
    context_summary = (
        f"Accessible public context matches: {context_audit['accessible_matches']}; compatible usable profiles: {context_audit['usable']}. Context assemblies were not requested."
        if context_audit
        else "Imported context rows are retained as supplied; this stage discovers no public context."
    )
    rows = []
    for row in audit["samples"]:
        calls = (
            "; ".join(
                f"{profile['called_loci']} / {profile['loci']} called ({profile['scheme_id']})"
                for profile in row["profiles"]
            )
            or "Unavailable"
        )
        reasons = "; ".join(row["reasons"]) or "Frozen profile and lineage evidence available"
        rows.append(
            "<tr>"
            + "".join(
                f"<td>{text(value)}</td>"
                for value in (
                    row["label"],
                    row["role"],
                    (f"ST{row['mlst_st']} ({row.get('mlst_scheme') or 'scheme unspecified'}; "
                     f"{row['mlst_status']})" if row['mlst_st_ready'] else 'Unresolved'),
                    calls,
                    "Available" if row["profile_ready"] else "Unavailable",
                    reasons,
                    ", ".join(row["missing_metadata"]) or "None",
                    ("Not bundled; retrieval not assessed"
                     if row["assembly_status"] == "retrieval_not_assessed"
                     else row["assembly_status"]),
                )
            )
            + "</tr>"
        )
    return (
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width, initial-scale=1">'
        "<title>ChronoClade preparation audit</title><style>"
        + report_styles()
        + '</style></head><body><main class="shell"><header class="identity">'
        '<div class="identity-mark">ChronoClade · PREPARATION AUDIT</div>'
        '<div class="identity-copy"><h1>Prepared dataset</h1>'
        "<p>Frozen evidence is saved in a validated portable bundle. Missing profiles, "
        "lineage assignments and metadata remain visible.</p></div></header>"
        '<section class="stage"><div class="stage-body"><h2>Prepared evidence</h2>'
        + parallel_summary
        +
        f"<p>Provider retrieval: {text(_operation_label(audit['operations']['provider_retrieval']))}. "
        f"{typing_label}: {text(_operation_label(audit['operations']['typing']))}. "
        f"cgMLST reference database: {text(_operation_label(audit['operations']['reference_database']))}.</p>"
        f'<div class="measure-strip"><div><span>Samples</span><b>{len(rows)}</b></div>'
        f"<div><span>cgMLST profiles available</span><b>{ready}</b></div>"
        f"<div><span>cgMLST profiles missing</span><b>{len(rows) - ready}</b></div></div>"
        f"<p>Missing metadata does not exclude a sample. {text(context_summary)}</p>"
        '<p>Optional lineage annotations are listed in the reasons below; missing '
        'cgLIN/HierCC labels do not make an available cgMLST profile unavailable. '
        'Public-context ST requirements and profile coverage are assessed separately.</p>'
        f"<h3>Can I continue?</h3><p>{text(eligibility)}</p>"
        '<div class="table-scroll" tabindex="0" role="region" aria-label="Sample readiness">'
        "<table><thead><tr>"
        + "".join(
            f'<th scope="col">{name}</th>'
            for name in (
                "Sample",
                "Role",
                "Classical MLST",
                "Profile calls",
                "Profile availability",
                "Reason",
                "Missing metadata",
                "Assembly",
            )
        )
        + "</tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></div>"
        '<details class="evidence-files" open><summary>Saved evidence</summary>'
        '<ul><li><a href="dataset.json">Dataset manifest</a></li>'
        '<li><a href="readiness.json">Machine-readable readiness audit</a></li>'
        '<li><a href="samples.csv">Samples and metadata</a></li>'
        '<li><a href="lineages.csv">Complete imported lineage evidence</a></li>'
        '<li><a href="provenance.csv">Field provenance</a></li></ul></details>'
        "</div></section></main></body></html>"
    )


def _next_stage_summary(audit):
    """Separate usable profile evidence from provider-specific typing completeness."""
    inputs = [row for row in audit["samples"] if row["role"] == "input"]
    missing = [row["sample_id"] for row in inputs if not row["profile_ready"]]
    incomplete = [row["sample_id"] for row in inputs
                  if any(not profile["catalogue_complete"] for profile in row["profiles"])]
    if missing:
        return (
            f'{len(missing)} of {len(inputs)} input genomes lack called cgMLST profiles. '
            'Supply compatible profiles or configure native typing and its reference database '
            'before interpreting a complete input-genome analysis. See the per-sample reasons below '
            '(also saved in readiness.json).'
        )
    text = (
        f'Called cgMLST profiles are available for all {len(inputs)} input genomes. '
        'Continue with cgmlst using the prepared dataset; compatible species and schemes are analysed separately. '
        'Missing cgLIN/HierCC labels do not prevent profile-distance analysis. '
    )
    if incomplete:
        text += (
            'The imported locus catalogue is incomplete: distances use the recorded loci only. '
            'A validated complete catalogue is required to interpret missing calls across the full scheme. '
        )
    text += (
        'Public-context discovery has separate provider-category and credential requirements. '
        'Missing dates can make time analysis unsupported; missing metadata does not remove genomes.'
    )
    return text


def run_prepare(
    input_path: str | Path,
    output: str | Path,
    *,
    catalogues: str | Path | None = None,
    input_kind: str = "auto",
    species: str | None = None,
    metadata: str | Path | None = None,
    typing_config: str | Path | None = None,
    query_typing: str | Path | None = None,
    public_typing: str | Path | None = None,
    cglin_export: str | Path | None = None,
    enrich_metadata: bool = True,
    client=None,
    _provider_sources: Path | None = None,
    _infer_mlst: bool = True,
    resume: bool = False,
) -> PrepareResult:
    """Import frozen canonical JSON or an existing bundle into a new stage output.

    Explicit catalogues are mandatory for canonical records, even when all calls
    are missing. Existing outputs are immutable. Local assembly references are
    copied and checksummed, without provider acquisition or context selection.
    """
    target = Path(output).expanduser().absolute()
    if target.exists():
        raise DatasetError(f"Preparation output already exists: {target}")
    if (
        input_kind == "auto"
        and isinstance(input_path, str)
        and ("://" in input_path or re.fullmatch(r"[1-9A-HJ-NP-Za-km-z]{20,22}", input_path))
    ):
        input_kind = "collection"
    if resume and input_kind not in {"collection", "accessions", "assemblies"}:
        raise DatasetError("--resume is supported only for live collection/accession/assembly preparation")
    if resume and not typing_config:
        raise DatasetError("--resume requires --typing-config so native typing work can be recovered")
    if input_kind in {"collection", "accessions", "assemblies"}:
        from dataclasses import asdict

        from chronoclade.prepare_provider import catalogues_for_records, resolve_queries
        from chronoclade.prepare_recovery import (
            DirectoryLock,
            RecoveryError,
            begin_job,
            bind_resolved_cohort,
            owned_directory,
            path_identity,
            save_job,
        )
        from chronoclade.execution_provenance import execution_provenance

        target.parent.mkdir(parents=True, exist_ok=True)
        work = target.with_name(target.name + ".work")
        input_file = Path(input_path).expanduser() if isinstance(input_path, (str, Path)) else None
        if input_file is not None and input_file.exists():
            input_identity = path_identity(input_file)
        else:
            input_identity = {"value": str(input_path)}
        binding = {
            "target": str(target),
            "input_kind": input_kind,
            "input": input_identity,
            "species": species,
            "sources": {
                "metadata": path_identity(metadata) if metadata else None,
                "typing_config": path_identity(typing_config) if typing_config else None,
                "query_typing": path_identity(query_typing) if query_typing else None,
                "public_typing": path_identity(public_typing) if public_typing else None,
                "cglin_export": path_identity(cglin_export) if cglin_export else None,
                "catalogues": path_identity(catalogues) if catalogues else None,
            },
            "flags": {"enrich_metadata": bool(enrich_metadata)},
            "provider_base_url": getattr(client, "base_url", None),
            "execution": execution_provenance(
                "native-prepare-recovery", "live-local-prepare-recovery-v1"
            ),
        }
        command = ["chronoclade", "prepare", str(input_path), "--out", str(target),
                   "--input-kind", input_kind]
        for option, value in (("--species", species), ("--metadata", metadata),
                              ("--typing-config", typing_config), ("--query-typing", query_typing),
                              ("--public-typing", public_typing), ("--cglin-export", cglin_export),
                              ("--catalogues", catalogues)):
            if value is not None:
                command.extend((option, str(value)))
        command.append("--enrich-metadata" if enrich_metadata else "--no-enrich-metadata")
        command.append("--resume")
        resume_command = shlex.join(command)
        lock_path = target.parent / f".{target.name}.prepare-recovery.lock"
        with DirectoryLock(lock_path):
            try:
                manifest = begin_job(work, binding, resume=resume)
            except RecoveryError as error:
                raise DatasetError(str(error)) from error
            attempt_id = uuid.uuid4().hex
            attempt = owned_directory(work, Path("attempts") / attempt_id, create=True)
            attempt_record = {
                "attempt_id": attempt_id,
                "status": "running",
                "evidence_directory": f"attempts/{attempt_id}",
            }
            manifest.setdefault("attempts", []).append(attempt_record)
            save_job(work, manifest)
            try:
                rows, provenance, conflicts, ledger, provider, operations = resolve_queries(
                    input_path,
                    attempt,
                    input_kind=input_kind,
                    species=species,
                    metadata=metadata,
                    client=client,
                    typing_config=typing_config,
                    query_typing=query_typing,
                    public_typing=public_typing,
                    cglin_export=cglin_export,
                    enrich_metadata=enrich_metadata,
                    native_checkpoint=work / "native_typing_cache",
                    native_attempt=owned_directory(attempt, "native_typing_attempts", create=True),
                    resume_typing=resume,
                    cohort_validator=lambda resolved: bind_resolved_cohort(work, manifest, resolved),
                )
                snapshot = provider.get("collection_snapshot")
                if snapshot:
                    snapshot_path = Path(snapshot).expanduser().resolve()
                    try:
                        snapshot_relative = snapshot_path.relative_to(attempt.resolve())
                    except ValueError as error:
                        raise DatasetError(
                            "Provider collection snapshot is outside the retained prepare sources"
                        ) from error
                    if not snapshot_path.is_file():
                        raise DatasetError("Provider collection snapshot is missing before publication")
                    provider["collection_snapshot"] = (
                        Path("sources") / "provider" / snapshot_relative
                    ).as_posix()
                explicit = _catalogues(Path(catalogues)) if catalogues else ()
                values = catalogues_for_records(rows, explicit)
                reference_parameters = {}
                if typing_config and operations.get("reference_database") == "ready_validated_configuration":
                    from chronoclade.pasteur_profile_normalization import reference_descriptor
                    reference_parameters["native_typing_reference"] = reference_descriptor(typing_config)
                source = attempt / "live_profiles.json"
                write_json(
                    source,
                    {
                        "records": rows,
                        "provenance": provenance,
                        "conflicts": conflicts,
                        "retrieval": ledger,
                        "parameters": {"provider_prepare": provider, "prepare_operations": operations,
                                       **reference_parameters},
                    },
                )
                # Preserve successful evidence from earlier failed attempts without
                # flattening paths or overwriting the current provider receipt.
                prior_attempts = [
                    item for item in manifest.get("attempts", [])
                    if item.get("attempt_id") != attempt_id
                ]
                if prior_attempts:
                    prior_root = attempt / "prior_attempts"
                    for prior in prior_attempts:
                        prior_source = work / prior["evidence_directory"]
                        prior_destination = prior_root / prior["attempt_id"]
                        for evidence in prior_source.rglob("*"):
                            if evidence.is_file() and evidence.suffix.lower() in {
                                ".json", ".jsonl", ".npy", ".csv", ".bin", ".raw", ".tsv",
                            }:
                                relative = evidence.relative_to(prior_source)
                                destination = prior_destination / relative
                                destination.parent.mkdir(parents=True, exist_ok=True)
                                shutil.copyfile(evidence, destination)
                if not values:
                    dataset = from_profile_records(
                        rows,
                        (),
                        provenance=provenance,
                        conflicts=conflicts,
                        retrieval=ledger,
                        parameters={"provider_prepare": provider, "prepare_operations": operations,
                                    **reference_parameters},
                    )
                    provisional = attempt / "sample_only"
                    manifest_path = write_dataset(dataset, provisional)
                    result = run_prepare(
                        manifest_path, target, input_kind="dataset", _provider_sources=attempt
                    )
                else:
                    catalogue_file = attempt / "live_catalogues.json"
                    write_json(catalogue_file, [asdict(c) for c in values])
                    result = run_prepare(
                        source, target, catalogues=catalogue_file, _provider_sources=attempt
                    )
                _validate_previous_stage(target)
                load_dataset(target / "dataset.json", mmap_mode=None)
                attempt_record["status"] = "published"
                manifest["status"] = "complete"
                save_job(work, manifest)
                # The final bundle has passed its sealed artifact validation. Remove
                # only this positively identified, job-owned scratch tree.
                saved = _json(work / "prepare-recovery.json")
                if (saved.get("binding") != binding or saved.get("status") != "complete"
                        or saved.get("job_id") != manifest.get("job_id")):
                    raise DatasetError("Recovery ownership marker changed; preserved work directory")
                shutil.rmtree(work)
                return result
            except KeyboardInterrupt as error:
                attempt_record["status"] = "interrupted"
                attempt_record["failure_type"] = "KeyboardInterrupt"
                manifest["status"] = "incomplete"
                if not typing_config:
                    shutil.rmtree(work, ignore_errors=True)
                    raise DatasetError("Preparation interrupted") from error
                save_job(work, manifest)
                raise DatasetError(
                    f"Preparation interrupted. Recovery checkpoint retained at {work}.\n"
                    f"Resume with: {resume_command}"
                ) from error
            except Exception as error:
                attempt_record["status"] = "failed"
                attempt_record["failure_type"] = type(error).__name__
                manifest["status"] = "incomplete"
                if not typing_config:
                    # Non-native live preparation keeps its historical ephemeral boundary.
                    shutil.rmtree(work, ignore_errors=True)
                    raise DatasetError(str(error)) from error
                save_job(work, manifest)
                detail = str(error) or f"{type(error).__name__}: preparation failed"
                raise DatasetError(
                    f"{detail}\nRecovery checkpoint retained at {work}.\nResume with: {resume_command}"
                ) from error
    if input_kind not in {"auto", "profiles", "dataset"}:
        raise DatasetError(
            "prepare requires profiles, dataset, collection, accessions or assemblies"
        )
    source = Path(input_path).expanduser().resolve()
    if species is not None and (not isinstance(species, str) or not species.strip()):
        raise DatasetError("--species must be nonempty text")
    if source.is_dir():
        source = source / "dataset.json"
    raw = _json(source)
    kind = input_kind
    if kind == "auto":
        kind = (
            "dataset" if isinstance(raw, dict) and raw.get("schema") == SCHEMA_NAME else "profiles"
        )
    catalogue_path = Path(catalogues).expanduser().resolve() if catalogues else None
    if kind == "dataset":
        if catalogue_path or species:
            raise DatasetError(
                "Existing dataset import preserves species/catalogues; omit overrides"
            )
        _validate_previous_stage(source.parent)
        dataset = load_dataset(source, mmap_mode=None)
    else:
        if catalogue_path is None:
            raise DatasetError(
                "Frozen profiles require --catalogues with an explicit locus catalogue"
            )
        dataset = _frozen_dataset(raw, _catalogues(catalogue_path), species)
    metadata_path = Path(metadata).expanduser().resolve() if metadata else None
    if metadata_path:
        dataset = _merge_frozen_metadata(dataset, metadata_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(prefix=f".{target.name}.prepare-", dir=target.parent))
    try:
        assembly_staging = scratch / "assembly_files"
        assembly_staging.mkdir()
        dataset, assemblies = _portable_assemblies(dataset, source.parent, assembly_staging)
        source_dataset_id = dataset.dataset_id
        from chronoclade.mlst_typing import fill_missing_mlst
        if _infer_mlst:
            typed_samples, mlst_audit = fill_missing_mlst(
                dataset.samples, assembly_staging, scratch / 'classical_mlst')
        else:
            typed_samples = dataset.samples
            mlst_audit = dataset.parameters.get('classical_mlst', {}).get('samples', [])
            if mlst_audit:
                write_json(scratch / 'classical_mlst' / 'audit.json', {
                    'schema': 'chronoclade.classical-mlst', 'schema_version': 1,
                    'policy': 'checked per-assembly results; collector inference disabled',
                    'samples': mlst_audit,
                })
        new_mlst_ids = ({row['sample_id'] for row in mlst_audit if row['status'] == 'resolved'}
                        if _infer_mlst else set())
        mlst_provenance = [
            {'sample_id': row['sample_id'], 'field': field, 'source': 'classical_mlst',
             'value': row[field], 'selected': True}
            for row in typed_samples if row['sample_id'] in new_mlst_ids
            for field in ('mlst_st', 'mlst_scheme', 'mlst_assembly_sha256', 'mlst_database_sha256')
        ]
        dataset = replace(dataset, samples=typed_samples,
                          provenance=(*dataset.provenance, *mlst_provenance),
                          retrieval=(*dataset.retrieval, *(
                              {'sample_id': row['sample_id'], 'operation': 'classical_mlst',
                               'status': row['status'], 'reason': row.get('reason')}
                              for row in mlst_audit if _infer_mlst and row['status'] in {'resolved', 'unresolved', 'failed'})),
                          dataset_id=None)
        source_evidence = {"path": "sources/input.json", "sha256": file_sha256(source)}
        parameters = dict(dataset.parameters)
        parameters['classical_mlst'] = {'policy': ('reuse supplied ST; exact primary-scheme calls for missing ST'
                                                 if _infer_mlst else 'checked per-assembly results; collector inference disabled'),
                                      'samples': mlst_audit,
                                      'evidence': 'sources/classical_mlst/audit.json' if mlst_audit else None}
        if typing_config:
            from chronoclade.pasteur_profile_normalization import reference_descriptor
            parameters["native_typing_reference"] = reference_descriptor(typing_config)
        parameters["prepare"] = {
            "input_kind": kind,
            "source": source_evidence,
            "source_dataset_id": source_dataset_id,
            "species": species,
            "assemblies": assemblies,
        }
        if catalogue_path:
            parameters["prepare"]["catalogues"] = {
                "path": "sources/catalogues.json",
                "sha256": file_sha256(catalogue_path),
            }
        if metadata_path:
            parameters["prepare"]["metadata"] = {
                "path": "sources/metadata.csv", "sha256": file_sha256(metadata_path),
            }
        dataset = replace(dataset, parameters=parameters, dataset_id=None)
        bundle = scratch / "bundle"
        manifest = write_dataset(dataset, bundle)
        sources = bundle / "sources"
        sources.mkdir()
        if (scratch / 'classical_mlst').exists():
            shutil.copytree(scratch / 'classical_mlst', sources / 'classical_mlst')
        shutil.copyfile(source, sources / "input.json")
        if catalogue_path:
            shutil.copyfile(catalogue_path, sources / "catalogues.json")
        if metadata_path:
            shutil.copyfile(metadata_path, sources / "metadata.csv")
        if (assembly_staging / "assemblies").exists():
            shutil.move(str(assembly_staging / "assemblies"), bundle / "assemblies")
        if _provider_sources:
            provider_sources = sources / "provider"
            provider_sources.mkdir()
            for evidence in _provider_sources.rglob("*"):
                if evidence.is_file() and evidence.suffix.lower() in {
                    ".json",
                    ".jsonl",
                    ".npy",
                    ".csv",
                    ".bin",
                    ".raw",
                    ".tsv",
                }:
                    if evidence.suffix.lower() == ".json":
                        _json(evidence)  # Reject credentials before snapshot publication.
                    if evidence.suffix.lower() == ".jsonl":
                        for line in evidence.read_text().splitlines():
                            check_json(json.loads(line), location=evidence.name)
                    relative = evidence.relative_to(_provider_sources)
                    destination = provider_sources / relative
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copyfile(evidence, destination)
        audit = _readiness(dataset, manifest)
        write_json(bundle / "readiness.json", audit)
        parallel_prepare = dataset.parameters.get('parallel_prepare')
        report = (_report(audit, parallel_prepare=parallel_prepare, evidence_root=bundle)
                  if parallel_prepare else _report(audit))
        (bundle / "report.html").write_text(report, encoding="utf-8")
        from chronoclade.stage_readme import write_stage_readme
        write_stage_readme(bundle, 'prepare', report='report.html',
                           readiness=audit['readiness'], published_directory=target,
                           readiness_summary=_next_stage_summary(audit))
        # Validate the actual saved dataset, then seal all stage-owned artifacts.
        load_dataset(manifest)
        artifacts = []
        for path in sorted(bundle.rglob("*")):
            if path.is_file():
                artifacts.append(
                    {
                        "path": path.relative_to(bundle).as_posix(),
                        "sha256": file_sha256(path),
                        "bytes": path.stat().st_size,
                    }
                )
        write_json(
            bundle / "prepare.json",
            {
                "schema": "chronoclade.prepare-result",
                "schema_version": 1,
                "software": {"name": "chronoclade", "version": __version__},
                "status": "complete",
                "input_kind": kind,
                "dataset_manifest": "dataset.json",
                "report": "report.html",
                "readiness_audit": "readiness.json",
                "readiness": audit["readiness"],
                "sample_ids": list(dataset.sample_ids),
                "artifacts": artifacts,
            },
        )
        bundle.rename(target)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    return PrepareResult(
        target / "dataset.json",
        target / "report.html",
        target / "readiness.json",
        target / "prepare.json",
    )
