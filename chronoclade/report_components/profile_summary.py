"""Compact dataset and subsample summaries for profile reports."""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from html import escape
from pathlib import Path, PurePosixPath
from typing import Any


def _mapping(value: object) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _items(value: object) -> list[Any]:
    return value if isinstance(value, list) else []


def _display(value: object) -> str:
    if value is None or value == "":
        return "Not reported"
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def _table(rows: list[dict[str, Any]], columns: list[tuple[str, str]], *, empty: str) -> str:
    if not rows:
        return f'<p class="empty">{escape(empty)}</p>'
    headings = "".join(f'<th scope="col">{escape(label)}</th>' for _, label in columns)
    body = "".join(
        "<tr>" + "".join(f"<td>{escape(_display(row.get(key)))}</td>" for key, _ in columns) + "</tr>"
        for row in rows
    )
    return (
        '<div class="table-wrap" tabindex="0" role="region" aria-label="Dataset summary table">'
        f"<table><thead><tr>{headings}</tr></thead><tbody>{body}</tbody></table></div>"
    )


def _origin(row: dict[str, Any]) -> str:
    value = str(row.get("origin", row.get("role", ""))).casefold()
    if value in {"local", "query", "focal", "input", "inputs"}:
        return "input"
    if value in {"context", "public", "comparison", "contextual"}:
        return "context"
    return "other"


def _known_count(counts: dict[str, int], field: str) -> int | str:
    return "Not fully classified" if counts["other"] else counts[field]


def _reported(value: object) -> bool:
    if value is None:
        return False
    text = str(value).strip()
    return text.casefold() not in {"", "unknown", "unavailable", "not reported", "none", "null", "n/a", "na", "—", "–", "-"}


def _geography_rows(data: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    aggregate: dict[tuple[str, tuple[str, ...]], dict[str, int]] = defaultdict(lambda: {"input": 0, "context": 0, "other": 0})
    countries: dict[str, dict[str, int]] = defaultdict(lambda: {"input": 0, "context": 0, "other": 0})
    for row in _items(data.get("metadata_geography")):
        if not isinstance(row, dict):
            continue
        count = row.get("count")
        if type(count) is not int or count < 0:
            continue
        country = str(row.get("country") or "Unknown")
        role = _origin(row)
        countries[country][role] += count
        region = row.get("region")
        nuts2 = row.get("nuts2")
        if _reported(region):
            aggregate[("region", (country, str(region).strip()))][role] += count
        if _reported(nuts2):
            location = (country, str(region).strip(), str(nuts2).strip()) if _reported(region) else (country, str(nuts2).strip())
            aggregate[("nuts2", location)][role] += count

    detail = []
    for (level, location), counts in sorted(aggregate.items()):
        detail.append({
            "level": {"region": "Region", "nuts2": "NUTS2"}[level],
            "place": " · ".join(location),
            "inputs": _known_count(counts, "input"),
            "context": _known_count(counts, "context"),
            "total": sum(counts.values()),
        })
    overall = [
        {"country": country, "inputs": _known_count(counts, "input"), "context": _known_count(counts, "context"),
         "total": sum(counts.values())}
        for country, counts in sorted(countries.items())
    ]
    return detail, overall


def _metadata_coverage(data: dict[str, Any]) -> str:
    summary = data.get("metadata_summary")
    if not summary:
        return ""
    rows: list[dict[str, Any]] = []
    if isinstance(summary, dict):
        for field, values in summary.items():
            if isinstance(values, dict):
                rows.append({"field": {"host": "Host", "collection_date": "Collection date"}.get(field, field), **values})
            else:
                rows.append({"field": field, "coverage": values})
    elif isinstance(summary, list):
        rows = [dict(row) for row in summary if isinstance(row, dict)]
    if not rows:
        return ""
    columns = [("field", "Metadata"), ("total", "Samples"), ("reported", "Reported"),
               ("missing", "Missing"), ("coverage", "Coverage")]
    present = [column for column in columns if column[0] == "field" or any(column[0] in row and _reported(row[column[0]]) for row in rows)]
    return '<h3>Host and collection-date coverage</h3>' + _table(
        rows, present, empty="No metadata coverage summary is available."
    )


def combined_dataset_summary(data: object) -> str:
    """Render geography and optional metadata coverage from a combined dataset."""
    payload = _mapping(data)
    geography, countries = _geography_rows(payload)
    output = ['<section class="profile-subsection" id="combined-dataset-summary">']
    if countries:
        output.append('<h3>Country totals</h3>')
        output.append(_table(
            countries,
            [("country", "Country"), ("inputs", "Inputs"), ("context", "Context"),
             ("total", "Total")],
            empty="No country totals are available.",
        ))
        detail = geography
        if detail:
            output.append('<details class="coverage-details"><summary>Region and NUTS2 detail</summary>')
            output.append(_table(
                detail,
                [("level", "Level"), ("place", "Reported location"), ("inputs", "Inputs"),
                 ("context", "Context"), ("total", "Total")],
                empty="No region or NUTS2 metadata is available.",
            ))
            output.append("</details>")
    else:
        output.append('<p class="empty">Geographic composition was not reported.</p>')
    for field in ("Host", "Collection year"):
        rows = [row for row in _items(payload.get("metadata_composition")) if isinstance(row, dict) and row.get("field") == field]
        if rows:
            output.append(f'<details class="coverage-details"><summary>{escape(field)} breakdown</summary>')
            output.append(_table(rows, [("value", field), ("input", "Inputs"), ("context", "Context"), ("total", "Total")], empty="No metadata breakdown is available."))
            output.append('</details>')
    output.append(_metadata_coverage(payload))
    output.append("</section>")
    return "".join(output)


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None
    return value if isinstance(value, dict) else None


def _selection_path(bundle: Path, run: dict[str, Any]) -> Path | None:
    ref = run.get("path")
    if not isinstance(ref, str) or not ref or "\\" in ref:
        return None
    relative = PurePosixPath(ref)
    if relative.is_absolute() or ".." in relative.parts:
        return None
    candidate = bundle / Path(*relative.parts)
    return candidate if candidate.is_file() else None


def _dataset_metadata(bundle: Path, ensemble: dict[str, Any], data: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Read sample metadata from the immutable prepared dataset when present."""
    direct = data.get("sample_metadata", data.get("samples"))
    if isinstance(direct, list):
        return {str(row["sample_id"]): row for row in direct
                if isinstance(row, dict) and row.get("sample_id")}
    source = _mapping(ensemble.get("source"))
    reference = _mapping(source.get("dataset"))
    raw_path = reference.get("path")
    if not isinstance(raw_path, str) or not raw_path:
        return {}
    dataset_path = (bundle / raw_path).resolve()
    if not dataset_path.is_file():
        return {}
    try:
        from chronoclade.datasets import load_dataset

        dataset = load_dataset(dataset_path)
    except (ImportError, OSError, ValueError, KeyError, TypeError):
        return {}
    return {str(row["sample_id"]): row for row in dataset.samples if row.get("sample_id")}


def _metadata_comparison(selection: dict[str, Any], pool_ids: set[str], metadata: dict[str, dict[str, Any]]) -> str:
    selected_ids = set(_items(selection.get("selected_sample_ids")))
    if not metadata or not selected_ids.intersection(metadata) or not pool_ids.intersection(metadata):
        return '<p class="empty">Country, collection-date and host comparison is unavailable because sample metadata is absent.</p>'
    rows = []
    for field, label in (("country", "Country"), ("collection_date", "Collection date"), ("host", "Host")):
        def known(identifier: str, field=field) -> bool:
            record = metadata.get(identifier, {})
            if field == "collection_date":
                return any(record.get(name) for name in ("collection_date", "date_start", "date_end"))
            return bool(record.get(field))

        selected_known = sum(known(identifier) for identifier in selected_ids)
        pool_known = sum(known(identifier) for identifier in pool_ids)
        rows.append({
            "field": label,
            "selected": f"{selected_known} / {len(selected_ids)} reported",
            "pool": f"{pool_known} / {len(pool_ids)} reported",
        })
    coverage_html = '<h4>Metadata coverage: selected set and full pool</h4>' + _table(
        rows, [("field", "Field"), ("selected", "Selected"), ("pool", "Full pool")],
        empty="Metadata comparison is unavailable.",
    )

    composition = []
    for field, label in (("country", "Country"), ("collection_date", "Collection year"), ("host", "Host")):
        def value(identifier, field=field):
            raw = metadata.get(identifier, {}).get(field)
            if not raw or str(raw).strip().casefold() in {"unknown", "not provided", "none"}:
                return "Unknown"
            if field == "collection_date":
                return str(raw)[:4] if str(raw)[:4].isdigit() else "Unknown"
            return str(raw).strip().capitalize() if field == "host" else str(raw)
        pool = Counter(value(identifier) for identifier in pool_ids)
        selected = Counter(value(identifier) for identifier in selected_ids)
        composition.extend({"field": label, "value": key, "selected": selected[key], "pool": count}
                           for key, count in sorted(pool.items()))
    return coverage_html + _table(composition, [("field", "Metadata"), ("value", "Value"),
                                                ("selected", "Selected genomes"), ("pool", "Full pool genomes")],
                                       empty="Metadata composition is unavailable.")


def _method_text(selection: dict[str, Any]) -> str:
    ranking = str(selection.get("nearest_ranking") or "")
    if ranking.startswith("raw_"):
        nearest = "Required nearest neighbours were ranked by raw differences when available, retaining boundary ties."
    elif ranking == "supplied_distance":
        nearest = "Required nearest neighbours were ranked using the supplied distance, retaining boundary ties."
    else:
        nearest = "The manifest records required nearest neighbours; its ranking metric is not specified."
    decisions = _items(selection.get("decisions"))
    reasons = {reason for row in decisions if isinstance(row, dict)
               for reason in _items(row.get("reasons")) if isinstance(reason, str)}
    if any(reason == "genetic_group_time_region_representative" for reason in reasons):
        optional = "Remaining representatives were sampled from recorded genetic-group, date, region and host strata; remaining slots used distance-based selection."
    else:
        optional = "Remaining context genomes were selected by the recorded distance-based procedure."
    if any("user_requested" in _items(row.get("reasons")) for row in decisions if isinstance(row, dict)):
        optional += " User-requested genomes were retained."
    return nearest + " " + optional


def subsample_summary(directory: str | Path, data: object) -> str:
    """Summarise a saved selection ensemble without expanding sample-level rows."""
    root = Path(directory)
    bundle = root / "selection_bundle"
    ensemble_path = bundle / "ensemble.json"
    output = ['<section class="profile-subsection" id="subsample-summary">']
    if not ensemble_path.is_file():
        output.append('<p class="empty">No subsamples were prepared for this dataset.</p></section>')
        return "".join(output)
    ensemble = _read_json(ensemble_path)
    if ensemble is None:
        output.append('<p class="empty">The saved subsample bundle could not be read.</p></section>')
        return "".join(output)
    runs = [run for run in _items(ensemble.get("selections")) if isinstance(run, dict)]
    selections: list[tuple[dict[str, Any], str]] = []
    for run in runs:
        path = _selection_path(bundle, run)
        selection = _read_json(path) if path else None
        if selection is not None:
            link = "selection_bundle/" + PurePosixPath(str(run["path"])).as_posix()
            selections.append((selection, link))
    if not selections:
        output.append('<p class="empty">The bundle contains no readable selection manifests.</p></section>')
        return "".join(output)

    first = selections[0][0]
    input_count = len(_items(first.get("query_ids")))
    context_count = len(_items(first.get("selected_context_ids")))
    pool_count = input_count + len(_items(first.get("available_context_ids")))
    output.append(f'<p><strong>{pool_count} genomes in the profile pool. The first representative set contains {input_count} input genomes + {context_count} context genomes = {input_count + context_count} genomes for SNP-tree analysis.</strong></p>')
    output.append('<p>These are saved input selections for the separate tree and time modules. Those downstream analyses have not been run by this cgMLST step. Alternative sets retain the inputs and required neighbours while varying the remaining representatives.</p>')
    reserve = _mapping(ensemble.get("parameters")).get("representative_reserve", 0)
    if reserve:
        output.append(f'<p>The requested reserve allows {int(reserve)} optional representatives in addition to mandatory nearest neighbours. The context budget expands when needed to preserve both; the available pool can limit the number selected.</p>')
    rows = []
    for selection, link in selections:
        query_ids = _items(selection.get("query_ids"))
        context_ids = _items(selection.get("selected_context_ids"))
        selected_ids = _items(selection.get("selected_sample_ids"))
        rows.append({
            "selection": selection.get("selection_id", "Selection"),
            "inputs": len(query_ids),
            "context": len(context_ids),
            "total": len(selected_ids) if selected_ids else len(query_ids) + len(context_ids),
            "nearest": len(set(_items(selection.get("required_nearest_context_ids")))),
            "manifest": f'<a href="{escape(link, quote=True)}">Download manifest</a>',
        })
    columns = [("selection", "Selection"), ("inputs", "Inputs"), ("context", "Context"),
               ("total", "Total selected"), ("nearest", "Mandatory nearest"), ("manifest", "Manifest")]
    # The manifest link is inserted as a fixed, path-validated anchor.
    headings = "".join(f'<th scope="col">{escape(label)}</th>' for _, label in columns)
    body = "".join(
        "<tr>" + "".join(f"<td>{escape(_display(row.get(key)))}</td>" for key, _ in columns[:-1])
        + f'<td>{row["manifest"]}</td></tr>' for row in rows
    )
    output.append(
        '<div class="table-wrap" tabindex="0" role="region" aria-label="Subsample selection table">'
        f"<table><thead><tr>{headings}</tr></thead><tbody>{body}</tbody></table></div>"
    )

    selected_sets = [set(_items(selection.get("selected_sample_ids"))) for selection, _ in selections]
    distinct_count = len({tuple(sorted(values)) for values in selected_sets})
    cumulative = set().union(*(set(_items(selection.get("selected_context_ids"))) for selection, _ in selections))
    pairwise = _items(ensemble.get("pairwise_overlap"))
    pairwise_values = [row.get("jaccard") for row in pairwise if isinstance(row, dict) and isinstance(row.get("jaccard"), (int, float))]
    overlap = f"{len(pairwise_values)} pairwise comparisons"
    if pairwise_values:
        overlap += f"; Jaccard overlap {min(pairwise_values):.2f}–{max(pairwise_values):.2f}"
    if len(selections) > 1:
        output.append(
            '<dl class="metrics">'
            f"<div><dt>Distinct alternatives</dt><dd>{distinct_count} / {len(selections)}</dd></div>"
            f"<div><dt>Cumulative context coverage</dt><dd>{len(cumulative)} unique genomes</dd></div>"
            f"<div><dt>Overlap</dt><dd>{escape(overlap)}</dd></div></dl>"
        )
    output.append(f"<p>{escape(_method_text(selections[0][0]))}</p>")
    metadata = _dataset_metadata(bundle, _mapping(ensemble), _mapping(data))
    pool_ids = set(_items(selections[0][0].get("query_ids"))) | set(_items(selections[0][0].get("available_context_ids")))
    for selection, _ in selections:
        output.append(f'<details class="coverage-details"><summary>{escape(_display(selection.get("selection_id")))} metadata comparison</summary>')
        output.append(_metadata_comparison(selection, pool_ids, metadata))
        output.append("</details>")
    output.append("</section>")
    return "".join(output)
