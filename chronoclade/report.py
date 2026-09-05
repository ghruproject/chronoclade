"""Human-readable temporal-signal reports."""

from __future__ import annotations

import csv
import base64
import math
import re
import zipfile
from html import escape
from io import StringIO
from pathlib import Path


RATE_DETAIL = re.compile(
    r"--rate:\s*([+-]?[0-9.]+(?:e[+-]?\d+)?)(?:\s*\+/-\s*([0-9.]+(?:e[+-]?\d+)?))?",
    re.IGNORECASE,
)
R_SQUARED_DETAIL = re.compile(r"--r\^2:\s*([+-]?[0-9.]+(?:e[+-]?\d+)?)", re.IGNORECASE)
CHI_SQUARED_DETAIL = re.compile(r"--chi\^2:\s*([+-]?[0-9.]+(?:e[+-]?\d+)?)", re.IGNORECASE)


def assess_temporal_signal(
    temporal: dict[str, object], *, p_value_threshold: float
) -> dict[str, object]:
    """Return a stable verdict and reader-facing explanation."""

    observed = temporal["observed"]
    if not isinstance(observed, dict):
        raise ValueError("temporal observed metrics must be a mapping")
    rate = float(observed["rate"])
    requested = int(temporal["requested_randomisations"])
    successful = int(temporal["successful_randomisations"])
    p_value = float(temporal["p_value_r_squared"])

    if requested == 0:
        return {
            "code": "not_assessed",
            "label": "Temporal signal not assessed",
            "reason": "No date randomisations were requested, so no dated tree was generated.",
            "supported": False,
        }
    if successful != requested:
        return {
            "code": "not_assessed",
            "label": "Temporal signal assessment incomplete",
            "reason": (
                f"Only {successful} of {requested} date randomisations completed successfully. "
                "The temporal-signal gate was not evaluated."
            ),
            "supported": False,
        }
    if rate <= 0:
        return {
            "code": "not_supported",
            "label": "Temporal signal not supported",
            "reason": (
                "The estimated evolutionary rate was not positive. Sampling dates do not "
                "support reliable time scaling for this lineage."
            ),
            "supported": False,
        }
    if p_value > p_value_threshold:
        return {
            "code": "not_supported",
            "label": "Temporal signal not supported",
            "reason": (
                f"The observed root-to-tip fit was not stronger than randomised dates at the "
                f"predefined threshold (p={p_value:.3g}; threshold={p_value_threshold:.3g})."
            ),
            "supported": False,
        }
    return {
        "code": "supported",
        "label": "Temporal signal supported",
        "reason": (
            "The lineage has a positive estimated rate and its observed root-to-tip fit was "
            "stronger than expected after randomly permuting sampling dates."
        ),
        "supported": True,
    }


def _metric(temporal: dict[str, object], name: str) -> float:
    observed = temporal["observed"]
    if not isinstance(observed, dict):
        raise ValueError("temporal observed metrics must be a mapping")
    return float(observed[name])


def _format_rate(value: float) -> str:
    return f"{value:.3g} substitutions/site/year"


def _embedded_font_css() -> str:
    """Embed the report display face so a saved report remains fully offline."""

    font = Path(__file__).parent / "data" / "fonts" / "archivo-latin-variable.woff2"
    if not font.is_file():
        return ""
    encoded = base64.b64encode(font.read_bytes()).decode("ascii")
    return (
        "@font-face{font-family:'Archivo';font-style:normal;font-weight:100 900;"
        "font-display:swap;src:url(data:font/woff2;base64," + encoded + ") format('woff2');}"
    )


def _pyplot():
    """Load matplotlib with a non-interactive backend only when figures are written."""

    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot

    return pyplot


def _plot_style(axis: object) -> None:
    axis.set_facecolor("#ffffff")
    axis.grid(axis="y", color="#d7d9df", linewidth=0.7, alpha=0.75)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["left"].set_color("#343842")
    axis.spines["bottom"].set_color("#343842")
    axis.tick_params(colors="#343842", labelsize=9)


def write_randomisation_csv(temporal: dict[str, object], output: Path) -> Path:
    """Write the complete observed and permuted clock results as tidy CSV."""

    raw_randomised = temporal.get("randomised", [])
    randomised = [value for value in raw_randomised if isinstance(value, dict)]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["analysis", "permutation", "rate", "r_squared", "exceeds_observed"],
        )
        writer.writeheader()
        observed_rate = _metric(temporal, "rate")
        observed_r_squared = _metric(temporal, "r_squared")
        writer.writerow(
            {
                "analysis": "observed",
                "permutation": "",
                "rate": observed_rate,
                "r_squared": observed_r_squared,
                "exceeds_observed": "",
            }
        )
        for index, value in enumerate(randomised, start=1):
            r_squared = float(value["r_squared"])
            writer.writerow(
                {
                    "analysis": "date_randomisation",
                    "permutation": index,
                    "rate": float(value["rate"]),
                    "r_squared": r_squared,
                    "exceeds_observed": str(r_squared >= observed_r_squared).lower(),
                }
            )
    return output


def write_randomisation_png(temporal: dict[str, object], output: Path) -> Path:
    """Write a publication-ready raster companion to the SVG randomisation plot."""

    raw_randomised = temporal.get("randomised", [])
    randomised = [value for value in raw_randomised if isinstance(value, dict)]
    pyplot = _pyplot()
    figure, axes = pyplot.subplots(1, 2, figsize=(10.6, 4.5), constrained_layout=True)
    panels = (
        (
            axes[0],
            [float(value["r_squared"]) for value in randomised],
            _metric(temporal, "r_squared"),
            "Root-to-tip fit",
            "R²",
        ),
        (
            axes[1],
            [float(value["rate"]) for value in randomised],
            _metric(temporal, "rate"),
            "Clock rate",
            "Substitutions/site/year",
        ),
    )
    for axis, values, observed, title, label in panels:
        _plot_style(axis)
        if values:
            bins = min(20, max(5, round(math.sqrt(len(values)))))
            axis.hist(values, bins=bins, color="#2855a6", alpha=0.78, edgecolor="white")
        axis.axvline(observed, color="#d63c2f", linewidth=2.4, label="Observed")
        axis.set_title(title, loc="left", fontsize=12, fontweight="bold", color="#17191f")
        axis.set_xlabel(label, color="#343842")
        axis.set_ylabel("Permutations", color="#343842")
        axis.legend(frameon=False, fontsize=9)
    figure.suptitle(
        f"Date randomisation · empirical p={float(temporal['p_value_r_squared']):.3g}",
        x=0.01,
        ha="left",
        fontsize=14,
        fontweight="bold",
        color="#17191f",
    )
    figure.savefig(
        output,
        dpi=220,
        facecolor="#ffffff",
        metadata={"Title": "Date randomisation", "Author": "ChronoClade"},
    )
    pyplot.close(figure)
    return output


def write_root_to_tip_png(source: Path, output: Path) -> Path | None:
    """Render TreeTime's root-to-tip table as a readable PNG scatter plot."""

    if not source.is_file():
        return None
    with source.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(line for line in handle if not line.startswith("#")))
    points: list[tuple[float, float, str]] = []
    for row in rows:
        normalized = {str(key).strip(): value for key, value in row.items()}
        try:
            points.append(
                (
                    float(str(normalized["date"]).strip()),
                    float(str(normalized["root-to-tip distance"]).strip()),
                    str(normalized.get("name", "")).strip(),
                )
            )
        except (KeyError, TypeError, ValueError):
            continue
    if not points:
        return None
    dates = [point[0] for point in points]
    distances = [point[1] for point in points]
    mean_x = sum(dates) / len(dates)
    mean_y = sum(distances) / len(distances)
    denominator = sum((value - mean_x) ** 2 for value in dates)
    slope = (
        sum(
            (x_value - mean_x) * (y_value - mean_y)
            for x_value, y_value in zip(dates, distances, strict=True)
        )
        / denominator
        if denominator
        else 0.0
    )
    intercept = mean_y - slope * mean_x
    pyplot = _pyplot()
    figure, axis = pyplot.subplots(figsize=(9.4, 5.4), constrained_layout=True)
    _plot_style(axis)
    axis.scatter(dates, distances, s=34, color="#2855a6", alpha=0.78, linewidth=0)
    lower, upper = min(dates), max(dates)
    axis.plot(
        [lower, upper],
        [intercept + slope * lower, intercept + slope * upper],
        color="#d63c2f",
        linewidth=2.2,
    )
    axis.set_title(
        "Genetic divergence through sampling time",
        loc="left",
        fontsize=15,
        fontweight="bold",
        color="#17191f",
    )
    axis.set_xlabel("Collection date", color="#343842")
    axis.set_ylabel("Root-to-tip distance (substitutions/site)", color="#343842")
    figure.savefig(
        output,
        dpi=220,
        facecolor="#ffffff",
        metadata={"Title": "Root-to-tip regression", "Author": "ChronoClade"},
    )
    pyplot.close(figure)
    return output


def _clock_confidence(directory: Path) -> dict[str, object]:
    """Read TreeTime clock and node-date confidence outputs when available."""

    clock_path = directory / "timetree" / "molecular_clock.txt"
    dates_path = directory / "timetree" / "dates.tsv"
    result: dict[str, object] = {
        "method": "TreeTime maximum likelihood",
        "interval_level": "90% max-posterior region",
    }
    if clock_path.is_file():
        text = clock_path.read_text(encoding="utf-8")
        rate = RATE_DETAIL.search(text)
        r_squared = R_SQUARED_DETAIL.search(text)
        chi_squared = CHI_SQUARED_DETAIL.search(text)
        if rate:
            result["rate"] = float(rate.group(1))
            result["rate_std_dev"] = float(rate.group(2)) if rate.group(2) else None
        if r_squared:
            result["r_squared"] = float(r_squared.group(1))
        if chi_squared:
            result["chi_squared"] = float(chi_squared.group(1))
    node_rows: list[dict[str, object]] = []
    if dates_path.is_file():
        with dates_path.open(encoding="utf-8") as handle:
            reader = csv.reader(
                (line for line in handle if not line.startswith("#")), delimiter="\t"
            )
            for row in reader:
                if len(row) < 5:
                    continue
                try:
                    numeric, lower, upper = map(float, row[2:5])
                except ValueError:
                    continue
                node_rows.append(
                    {
                        "node": row[0],
                        "date": row[1],
                        "numeric_date": numeric,
                        "lower_bound": lower,
                        "upper_bound": upper,
                        "interval_width_years": upper - lower,
                    }
                )
        internal = [row for row in node_rows if str(row["node"]).startswith("NODE_")]
        if internal:
            root = internal[0]
            widths = sorted(float(row["interval_width_years"]) for row in internal)
            result.update(
                {
                    "root_node": root["node"],
                    "root_date": root["date"],
                    "root_numeric_date": root["numeric_date"],
                    "root_lower_bound": root["lower_bound"],
                    "root_upper_bound": root["upper_bound"],
                    "internal_nodes": len(internal),
                    "median_internal_interval_width_years": widths[len(widths) // 2],
                    "maximum_internal_interval_width_years": max(widths),
                }
            )
    result["node_rows"] = node_rows
    return result


def write_timetree_confidence(directory: Path) -> dict[str, object]:
    """Export TreeTime confidence summaries and per-node dates as CSV."""

    confidence = _clock_confidence(directory)
    node_rows = confidence.pop("node_rows", [])
    summary_path = directory / "timetree_confidence.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        for key, value in confidence.items():
            writer.writerow([key, "" if value is None else value])
    if node_rows:
        node_path = directory / "node_dates.csv"
        with node_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(node_rows[0]))
            writer.writeheader()
            writer.writerows(node_rows)
    return confidence


def write_timetree_figures(directory: Path, confidence: dict[str, object]) -> Path | None:
    """Render a calendar-time phylogeny with visible node-date intervals."""

    source = directory / "timetree" / "timetree.nexus"
    if not source.is_file():
        return None
    from Bio import Phylo

    node_dates = directory / "node_dates.csv"
    if not node_dates.is_file():
        return None
    with node_dates.open(newline="", encoding="utf-8") as handle:
        date_rows = {
            row["node"]: {
                "date": float(row["numeric_date"]),
                "lower": float(row["lower_bound"]),
                "upper": float(row["upper_bound"]),
            }
            for row in csv.DictReader(handle)
        }

    nexus_text = source.read_text(encoding="utf-8")
    tree_line = next(
        line for line in nexus_text.splitlines() if line.lstrip().casefold().startswith("tree ")
    )
    tree = Phylo.read(StringIO(tree_line.split("=", 1)[1].strip()), "newick")
    pyplot = _pyplot()
    tip_count = tree.count_terminals()
    figure, axis = pyplot.subplots(
        figsize=(12.2, max(7.0, min(14.0, tip_count * 0.15))), constrained_layout=True
    )
    terminals = tree.get_terminals()
    y_positions = {clade: float(index) for index, clade in enumerate(reversed(terminals))}

    def y_position(clade: object) -> float:
        if clade in y_positions:
            return y_positions[clade]
        children = list(clade.clades)
        value = sum(y_position(child) for child in children) / len(children)
        y_positions[clade] = value
        return value

    y_position(tree.root)

    def numeric_date(clade: object, parent_date: float | None = None) -> float:
        row = date_rows.get(str(clade.name))
        if row:
            return row["date"]
        if parent_date is None:
            return float(confidence.get("root_numeric_date", 0.0))
        return parent_date + float(clade.branch_length or 0.0)

    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal() or not clade.name or str(clade.name) not in date_rows:
            continue
        row = date_rows[str(clade.name)]
        axis.hlines(
            y_positions[clade],
            row["lower"],
            row["upper"],
            color="#e2aaa3",
            linewidth=4.2,
            alpha=0.72,
            zorder=1,
        )

    def draw_clade(clade: object, parent_date: float | None = None) -> None:
        x_value = numeric_date(clade, parent_date)
        children = list(clade.clades)
        if children:
            child_y = [y_positions[child] for child in children]
            axis.vlines(
                x_value,
                min(child_y),
                max(child_y),
                color="#17191f",
                linewidth=0.65,
                zorder=2,
            )
        for child in children:
            child_x = numeric_date(child, x_value)
            axis.hlines(
                y_positions[child],
                x_value,
                child_x,
                color="#17191f",
                linewidth=0.7,
                zorder=2,
            )
            draw_clade(child, x_value)

    draw_clade(tree.root)
    tip_x = [numeric_date(clade) for clade in terminals]
    tip_y = [y_positions[clade] for clade in terminals]
    axis.scatter(tip_x, tip_y, s=11, color="#2855a6", zorder=3, label="Sampled genome")
    root_x = numeric_date(tree.root)
    axis.scatter([root_x], [y_positions[tree.root]], s=32, color="#d63c2f", zorder=4, label="Root")
    root_lower = confidence.get("root_lower_bound")
    root_upper = confidence.get("root_upper_bound")
    if root_lower is not None and root_upper is not None:
        axis.annotate(
            f"Root {root_x:.1f}\n90% interval {float(root_lower):.1f}–{float(root_upper):.1f}",
            xy=(root_x, y_positions[tree.root]),
            xytext=(12, 12),
            textcoords="offset points",
            fontsize=9,
            color="#8f241d",
            fontweight="bold",
        )
    from matplotlib.lines import Line2D

    legend_handles = [
        Line2D([0], [0], color="#17191f", linewidth=1, label="Time-scaled branch"),
        Line2D([0], [0], color="#e2aaa3", linewidth=4.2, label="90% node-date interval"),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#2855a6",
            markeredgecolor="none",
            label="Sampled genome",
        ),
    ]
    axis.legend(handles=legend_handles, loc="upper left", frameon=False, ncol=3, fontsize=9)
    axis.set_xlabel("Calendar year")
    axis.set_yticks([])
    axis.set_ylabel("")
    axis.set_title(
        "Dated phylogeny with node-time uncertainty",
        loc="left",
        fontsize=15,
        fontweight="bold",
        pad=26,
    )
    axis.grid(axis="x", color="#d7d9df", linewidth=0.7, alpha=0.8)
    axis.margins(x=0.03, y=0.015)
    axis.spines["left"].set_visible(False)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.spines["bottom"].set_color("#343842")
    output = directory / "timetree_with_confidence.svg"
    figure.savefig(output, facecolor="white")
    figure.savefig(
        directory / "timetree.png",
        dpi=220,
        facecolor="white",
        metadata={
            "Title": "Dated phylogeny with node-time uncertainty",
            "Author": "ChronoClade",
        },
    )
    pyplot.close(figure)
    return output


def _histogram_panel(
    values: list[float],
    *,
    observed: float,
    x: float,
    title: str,
    x_label: str,
    formatter: str,
) -> str:
    panel_width = 430.0
    plot_x = x + 52
    plot_y = 72.0
    plot_width = 340.0
    plot_height = 220.0
    title_text = escape(title)
    label_text = escape(x_label)

    if not values:
        return f"""
        <g>
          <text x="{x + panel_width / 2}" y="36" text-anchor="middle" class="title">{title_text}</text>
          <rect x="{plot_x}" y="{plot_y}" width="{plot_width}" height="{plot_height}" class="plot-bg"/>
          <text x="{x + panel_width / 2}" y="185" text-anchor="middle" class="empty">No completed permutations</text>
          <text x="{x + panel_width / 2}" y="330" text-anchor="middle" class="axis-label">{label_text}</text>
        </g>
        """

    lower = min([*values, observed])
    upper = max([*values, observed])
    if math.isclose(lower, upper):
        padding = abs(lower) * 0.1 or 1.0
        lower -= padding
        upper += padding
    else:
        padding = (upper - lower) * 0.05
        lower -= padding
        upper += padding

    bin_count = min(20, max(5, round(math.sqrt(len(values)))))
    bin_width = (upper - lower) / bin_count
    counts = [0] * bin_count
    for value in values:
        index = min(bin_count - 1, max(0, int((value - lower) / bin_width)))
        counts[index] += 1
    maximum = max(counts) or 1
    bar_width = plot_width / bin_count
    bars: list[str] = []
    for index, count in enumerate(counts):
        height = plot_height * count / maximum
        bars.append(
            f'<rect x="{plot_x + index * bar_width + 1:.2f}" '
            f'y="{plot_y + plot_height - height:.2f}" width="{max(1, bar_width - 2):.2f}" '
            f'height="{height:.2f}" class="bar"/>'
        )

    observed_x = plot_x + ((observed - lower) / (upper - lower)) * plot_width
    lower_text = format(lower, formatter)
    upper_text = format(upper, formatter)
    observed_text = format(observed, formatter)
    return f"""
    <g>
      <text x="{x + panel_width / 2}" y="36" text-anchor="middle" class="title">{title_text}</text>
      <rect x="{plot_x}" y="{plot_y}" width="{plot_width}" height="{plot_height}" class="plot-bg"/>
      {"".join(bars)}
      <line x1="{observed_x:.2f}" y1="{plot_y - 8}" x2="{observed_x:.2f}" y2="{plot_y + plot_height}" class="observed"/>
      <text x="{observed_x:.2f}" y="{plot_y - 16}" text-anchor="middle" class="observed-label">Observed {escape(observed_text)}</text>
      <line x1="{plot_x}" y1="{plot_y + plot_height}" x2="{plot_x + plot_width}" y2="{plot_y + plot_height}" class="axis"/>
      <text x="{plot_x}" y="{plot_y + plot_height + 20}" text-anchor="start" class="tick">{escape(lower_text)}</text>
      <text x="{plot_x + plot_width}" y="{plot_y + plot_height + 20}" text-anchor="end" class="tick">{escape(upper_text)}</text>
      <text x="{x + panel_width / 2}" y="330" text-anchor="middle" class="axis-label">{label_text}</text>
    </g>
    """


def write_randomisation_plot(temporal: dict[str, object], output: Path) -> Path:
    """Write an accessible two-panel SVG comparing observed and permuted metrics."""

    raw_randomised = temporal.get("randomised", [])
    randomised = [value for value in raw_randomised if isinstance(value, dict)]
    r_squared = [float(value["r_squared"]) for value in randomised]
    rates = [float(value["rate"]) for value in randomised]
    observed_r_squared = _metric(temporal, "r_squared")
    observed_rate = _metric(temporal, "rate")
    p_value = float(temporal["p_value_r_squared"])
    requested = int(temporal["requested_randomisations"])
    successful = int(temporal["successful_randomisations"])

    left = _histogram_panel(
        r_squared,
        observed=observed_r_squared,
        x=20,
        title="Root-to-tip fit under randomised dates",
        x_label="R²",
        formatter=".3g",
    )
    right = _histogram_panel(
        rates,
        observed=observed_rate,
        x=470,
        title="Clock rate under randomised dates",
        x_label="Substitutions/site/year",
        formatter=".2e",
    )
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="920" height="390" viewBox="0 0 920 390" role="img" aria-labelledby="title description">
  <title id="title">Date-randomisation comparison</title>
  <desc id="description">Histograms compare the observed root-to-tip R squared and clock rate with values from permuted sampling dates.</desc>
  <style>
    text {{ font-family: Archivo, sans-serif; fill: #243340; }}
    .title {{ font-size: 16px; font-weight: 700; }}
    .axis-label {{ font-size: 13px; font-weight: 600; }}
    .tick, .empty {{ font-size: 12px; fill: #607080; }}
    .plot-bg {{ fill: #f5f8fa; stroke: #d6e0e6; }}
    .bar {{ fill: #6285a3; opacity: 0.82; }}
    .observed {{ stroke: #c6403d; stroke-width: 3; }}
    .observed-label {{ font-size: 12px; font-weight: 700; fill: #a52f2c; }}
    .axis {{ stroke: #607080; stroke-width: 1; }}
    .footer {{ font-size: 13px; fill: #40515e; }}
  </style>
  {left}
  {right}
  <text x="460" y="372" text-anchor="middle" class="footer">{successful} of {requested} permutations completed · root-to-tip randomisation p={p_value:.3g}</text>
</svg>
"""
    output.write_text(svg, encoding="utf-8")
    return output


def _file_link(relative_path: str, label: str) -> str:
    return f'<a href="{escape(relative_path, quote=True)}">{escape(label)}</a>'


def _download_list(items: list[tuple[str, str, str]], directory: Path) -> str:
    rows = []
    for relative_path, label, description in items:
        if not (directory / relative_path).is_file():
            continue
        suffix = Path(relative_path).suffix.removeprefix(".").upper() or "FILE"
        rows.append(
            '<li class="download"><a href="'
            f'{escape(relative_path, quote=True)}"><span>{escape(label)}</span>'
            f"<small>{escape(description)}</small><b>{escape(suffix)}</b></a></li>"
        )
    return f'<ul class="downloads">{"".join(rows)}</ul>' if rows else ""


def _table_region(table: str, *, label: str, compact: bool = False) -> str:
    """Wrap a wide result table in a labelled, keyboard-scrollable region."""

    classes = "table-scroll compact-table" if compact else "table-scroll"
    return (
        f'<div class="{classes}" tabindex="0" role="region" '
        f'aria-label="{escape(label, quote=True)}">'
        '<span class="scroll-hint">Scroll horizontally to inspect every column</span>'
        f"{table}</div>"
    )


def write_supporting_bundle(directory: Path) -> Path:
    """Package reader-facing evidence files without duplicating large intermediates."""

    bundle = directory / "supporting_results.zip"
    allowed = {".csv", ".tsv", ".json", ".txt", ".svg", ".png", ".nexus", ".newick"}
    excluded_parts = {"logs", "randomisations"}
    with zipfile.ZipFile(bundle, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(directory.rglob("*")):
            if (
                not path.is_file()
                or path == bundle
                or path.suffix.lower() not in allowed
                or excluded_parts.intersection(path.relative_to(directory).parts)
            ):
                continue
            archive.write(path, path.relative_to(directory))
    return bundle


def _display_summary(value: object) -> str:
    if not isinstance(value, dict) or not value.get("comparisons"):
        return "No comparisons"
    return (
        f"n={int(value['comparisons'])}; median {float(value['median']):.3g}; "
        f"range {int(value['minimum'])}–{int(value['maximum'])} SNPs"
    )


def _public_health_visual(report: dict[str, object], directory: Path) -> tuple[str, str]:
    """Return the reader summary and the detailed public-health evidence separately."""
    raw_public_health = report.get("public_health", {})
    public_health = raw_public_health if isinstance(raw_public_health, dict) else {}
    raw_scenario = public_health.get("scenario", {})
    scenario = raw_scenario if isinstance(raw_scenario, dict) else {}
    code = str(scenario.get("code", "indeterminate"))
    label = str(scenario.get("label", "Indeterminate"))
    confidence = str(scenario.get("confidence", "low"))
    reasons = scenario.get("reasons", [])
    actions = scenario.get("recommended_follow_up", [])
    evidence = scenario.get("evidence", [])
    if not isinstance(reasons, list):
        reasons = []
    if not isinstance(actions, list):
        actions = []
    if not isinstance(evidence, list):
        evidence = []
    evidence_rows = []
    for item in evidence:
        if not isinstance(item, dict):
            continue
        evidence_rows.append(
            "<tr>"
            f"<td>{escape(str(item.get('id', '')).replace('_', ' ').title())}</td>"
            f'<td><span class="evidence {escape(str(item.get("status", "unavailable")))}">'
            f"{escape(str(item.get('status', 'unavailable')).title())}</span></td>"
            f"<td>{escape(str(item.get('finding', '')))}</td>"
            "</tr>"
        )

    raw_distances = public_health.get("distance_summary", {})
    distances = raw_distances if isinstance(raw_distances, dict) else {}
    raw_categories = distances.get("categories", {})
    categories = raw_categories if isinstance(raw_categories, dict) else {}
    category_labels = (
        ("focal_focal", "Focal–focal"),
        ("focal_context", "Focal–context"),
        ("same_month", "Focal pairs in the same month"),
        ("same_year", "Focal pairs in the same year"),
        ("between_years", "Focal pairs between years"),
        ("within_candidate_groups", "Within candidate local groups"),
        ("between_candidate_groups", "Between candidate local groups"),
        ("same_patient", "Within-patient focal pairs"),
        ("different_patients", "Between-patient focal pairs"),
    )
    distance_rows = "".join(
        f"<tr><td>{escape(label_text)}</td><td>{escape(_display_summary(categories.get(key)))}</td></tr>"
        for key, label_text in category_labels
    )

    raw_topology = public_health.get("topology", {})
    topology = raw_topology if isinstance(raw_topology, dict) else {}
    raw_groups = topology.get("groups", [])
    groups = raw_groups if isinstance(raw_groups, list) else []
    group_rows = []
    for group in groups:
        if not isinstance(group, dict):
            continue
        sample_ids = group.get("sample_ids", [])
        if not isinstance(sample_ids, list):
            sample_ids = []
        locations = group.get("locations", [])
        if not isinstance(locations, list):
            locations = []
        nearest = group.get("nearest_context")
        nearest_text = "—"
        if isinstance(nearest, dict):
            nearest_text = (
                f"{nearest.get('sample_id', '')} ({nearest.get('clonal_snps', '—')} SNPs)"
            )
        group_rows.append(
            "<tr>"
            f"<td>{escape(str(group.get('group_id', '')))}</td>"
            f"<td>{int(group.get('sample_count', 0))}</td>"
            f"<td>{int(group.get('patient_count', 0)) or '—'}</td>"
            f"<td>{escape(str(group.get('first_collection_date', '') or '—'))} – "
            f"{escape(str(group.get('last_collection_date', '') or '—'))}</td>"
            f"<td>{escape(', '.join(str(value) for value in locations) or '—')}</td>"
            f"<td>{escape(_display_summary(group.get('within_group_clonal_snps')))}</td>"
            f"<td>{escape(nearest_text)}</td>"
            f"<td><details><summary>View IDs</summary>{escape(', '.join(str(value) for value in sample_ids))}</details></td>"
            "</tr>"
        )
    group_table = (
        _table_region(
            "<table><thead><tr><th>Candidate group</th><th>Focal genomes</th><th>Patients</th>"
            "<th>Date span</th><th>Locations</th><th>Within-group clonal distance</th>"
            "<th>Nearest final context</th><th>Samples</th></tr></thead>"
            f"<tbody>{''.join(group_rows)}</tbody></table>",
            label="Candidate local groups",
        )
        if group_rows
        else "<p>No topology-defined focal groups were available.</p>"
    )

    raw_sensitivity = public_health.get("patient_sensitivity", {})
    sensitivity = raw_sensitivity if isinstance(raw_sensitivity, dict) else {}
    excluded = sensitivity.get("excluded_repeated_patient_samples", [])
    if not isinstance(excluded, list):
        excluded = []
    sensitivity_text = str(
        sensitivity.get("interpretation", "Patient-level sensitivity was not available.")
    )
    if not bool(sensitivity.get("patient_metadata_complete")):
        sensitivity_text += " Patient identifiers are incomplete for the focal samples."
    if not excluded:
        exclusion_text = (
            "No repeated-patient samples would be excluded in the deterministic "
            "one-isolate-per-patient view."
        )
    elif len(excluded) == 1:
        exclusion_text = (
            "One repeated-patient sample would be excluded in the deterministic "
            "one-isolate-per-patient view."
        )
    else:
        exclusion_text = (
            f"{len(excluded)} repeated-patient samples would be excluded in the deterministic "
            "one-isolate-per-patient view."
        )

    heatmap = directory / "clonal_snp_heatmap.svg"
    heatmap_visual = (
        '<a class="figure-expand" href="clonal_snp_heatmap.svg" title="Open the full-resolution heatmap">'
        '<img src="clonal_snp_heatmap.svg" alt="Heatmap of recombination-filtered pairwise SNP counts"></a>'
        '<p><a href="clonal_snp_heatmap.svg">Open the full-resolution heatmap</a></p>'
        if heatmap.is_file()
        else '<p class="missing">The clonal SNP heatmap was not available.</p>'
    )
    evidence_table = (
        _table_region(
            "<table><thead><tr><th>Evidence</th><th>Status</th><th>Finding</th></tr></thead>"
            f"<tbody>{''.join(evidence_rows)}</tbody></table>",
            label="Evidence ledger",
        )
        if evidence_rows
        else "<p>No structured evidence ledger was available.</p>"
    )
    summary = f"""
  <section class="scenario {escape(code)}">
    <h3>What does the genomic evidence support?</h3>
    <strong>{escape(label)}</strong>
    <p><span class="confidence">{escape(confidence)} confidence</span></p>
    <h3>Why?</h3>
    <ul>{"".join(f"<li>{escape(str(reason))}</li>" for reason in reasons)}</ul>
    <p class="guardrail">{escape(str(scenario.get("guardrail", "This analysis does not establish direct transmission.")))}</p>
  </section>

  <section class="card action">
    <h3>What should happen next?</h3>
    <h4>Recommended follow-up</h4>
    <ul>{"".join(f"<li>{escape(str(action))}</li>" for action in actions)}</ul>
  </section>
"""

    technical = f"""
  <section class="card rule">
    <h3>Interpretation rule</h3>
    <p>{escape(str(scenario.get("decision_rule", "No scenario rule was available.")))}</p>
    <p>This automated summary is provisional and requires epidemiological review.</p>
  </section>

  <section class="card">
    <h3>Evidence ledger</h3>
    <p>Every part of the working interpretation is exposed here rather than combined into an opaque score.</p>
    {evidence_table}
  </section>

  <section class="card">
    <h3>Candidate local groups</h3>
    <p>{escape(str(topology.get("interpretation", "No topology summary was available.")))}</p>
    {group_table}
  </section>

  <section class="card">
    <h3>Recombination-filtered genomic distances</h3>
    <p>Counts use only pair-specific A/C/G/T sites after ClonalFrameML recombination filtering. No universal SNP threshold has been applied.</p>
    {heatmap_visual}
    {_table_region('<table><thead><tr><th>Comparison</th><th>Clonal SNP summary</th></tr></thead><tbody>' + distance_rows + '</tbody></table>', label='Recombination-filtered genomic distance summaries', compact=True)}
    <p>{_file_link("clonal_pairwise_distances.tsv", "Download pairwise SNPs and callable sites")} · {_file_link("clonal_snp_matrix.tsv", "Download SNP matrix")} · {_file_link("pairwise_callable_sites.tsv", "Download callable-site matrix")}</p>
  </section>

  <section class="card">
    <h3>Patient-level sensitivity</h3>
    <p>{escape(sensitivity_text)}</p>
    <p>{escape(exclusion_text)}</p>
  </section>
"""
    return summary, technical


def write_lineage_report(
    report: dict[str, object], *, directory: Path, p_value_threshold: float
) -> Path:
    """Write a standalone, reader-facing HTML report for one lineage."""

    temporal = report["temporal_signal"]
    if not isinstance(temporal, dict):
        raise ValueError("lineage report lacks temporal metrics")
    assessment = assess_temporal_signal(temporal, p_value_threshold=p_value_threshold)
    status = str(assessment["code"])
    raw_context = report.get("context", {})
    context = raw_context if isinstance(raw_context, dict) else {}
    context_count = int(context.get("context_samples", 0))
    context_locations = context.get("context_locations", [])
    if not isinstance(context_locations, list):
        context_locations = []
    public_health_summary, public_health_evidence = _public_health_visual(report, directory)
    observed_rate = _metric(temporal, "rate")
    observed_r_squared = _metric(temporal, "r_squared")
    p_value = float(temporal["p_value_r_squared"])
    successful = int(temporal["successful_randomisations"])

    root_plot = directory / "clock" / "root_to_tip_regression.svg"
    timetree_plot = directory / "timetree" / "timetree.svg"
    randomisation_plot = directory / "date_randomisation.svg"
    write_randomisation_plot(temporal, randomisation_plot)
    write_randomisation_png(temporal, directory / "date_randomisation.png")
    write_randomisation_csv(temporal, directory / "date_randomisation.csv")
    write_root_to_tip_png(directory / "clock" / "rtt.csv", directory / "root_to_tip.png")
    confidence = write_timetree_confidence(directory) if bool(assessment["supported"]) else {}
    confidence_plot: Path | None = None
    if confidence:
        confidence_plot = write_timetree_figures(directory, confidence)

    root_visual = (
        '<img src="clock/root_to_tip_regression.svg" alt="TreeTime root-to-tip regression plot">'
        if root_plot.is_file()
        else '<p class="missing">Root-to-tip plot was not available. See the clock output files.</p>'
    )
    timetree_visual = ""
    if bool(assessment["supported"]) and timetree_plot.is_file():
        displayed_timetree = (
            "timetree_with_confidence.svg"
            if confidence_plot is not None
            else "timetree/timetree.svg"
        )
        root_numeric = confidence.get("root_numeric_date")
        root_date = (
            f"{float(root_numeric):.1f}" if root_numeric is not None else "not available"
        )
        lower = confidence.get("root_lower_bound")
        upper = confidence.get("root_upper_bound")
        root_interval = (
            f"{float(lower):.1f}–{float(upper):.1f}"
            if lower is not None and upper is not None
            else "not available"
        )
        rate_std = confidence.get("rate_std_dev")
        rate_std_text = f" ± {float(rate_std):.2g}" if rate_std is not None else ""
        confidence_downloads = _download_list(
            [
                (
                    "timetree_with_confidence.svg",
                    "Publication dated phylogeny — SVG",
                    "Scalable figure with node-date intervals",
                ),
                (
                    "timetree.png",
                    "Publication dated phylogeny — PNG",
                    "High-resolution figure with node-date intervals",
                ),
                (
                    "timetree/timetree.svg",
                    "Raw TreeTime figure — SVG",
                    "Advanced: unmodified TreeTime output",
                ),
                ("timetree/timetree.nexus", "Dated tree", "Tree with time-scaled branch lengths"),
                (
                    "timetree_confidence.csv",
                    "Confidence summary",
                    "Method, rate uncertainty and interval widths",
                ),
                ("node_dates.csv", "Node dates", "All inferred node dates and 90% intervals"),
                ("timetree/dates.tsv", "TreeTime node-date output", "Original TreeTime table"),
                (
                    "timetree/molecular_clock.txt",
                    "TreeTime clock output",
                    "Original fitted clock statistics",
                ),
            ],
            directory,
        )
        timetree_visual = f"""
        <section class="stage" id="time-tree" data-stage="3">
          <div class="stage-index"><span>3</span><b>Time scale</b></div>
          <div class="stage-body">
            <div class="stage-head">
              <div><h2>Estimate the dated phylogeny</h2><p class="question">Only reached because the temporal-signal gate passed.</p></div>
              <strong class="decision proceed">PROCEED</strong>
            </div>
            <p>TreeTime now estimates calendar dates for internal nodes. The pale vermillion bars
            on the phylogeny show the 90% node-date interval for each inferred node. TreeTime
            refits the final dated tree, so its R² can differ slightly from the exploratory fit.</p>
            <div class="confidence-grid">
              <div><small>Clock rate</small><b>{float(confidence.get("rate", observed_rate)):.3g}{escape(rate_std_text)}</b><span>subs/site/year ± 1 SD</span></div>
              <div><small>Final dated-tree fit</small><b>R² {float(confidence.get("r_squared", observed_r_squared)):.3g}</b><span>after time-tree refitting</span></div>
              <div><small>Root estimate</small><b>{escape(str(root_date))}</b><span>calendar year; 90% interval {escape(root_interval)}</span></div>
              <div><small>Internal nodes</small><b>{int(confidence.get("internal_nodes", 0))}</b><span>with inferred dates</span></div>
              <div><small>Median interval width</small><b>{float(confidence.get("median_internal_interval_width_years", 0.0)):.2f} years</b><span>across internal nodes</span></div>
              <div><small>Widest interval</small><b>{float(confidence.get("maximum_internal_interval_width_years", 0.0)):.2f} years</b><span>least precise node date</span></div>
            </div>
            <div class="evidence-layout"><figure><div class="figure-scroll" tabindex="0" role="region" aria-label="Dated phylogeny; scroll horizontally to inspect detail"><a class="figure-expand" href="{displayed_timetree}" title="Open the full-resolution dated phylogeny"><img src="{displayed_timetree}" alt="Time-scaled phylogeny with visible 90% node-date intervals"></a></div><figcaption>Calendar-time phylogeny. Branches are black, sampled genomes are blue, and pale vermillion bars are 90% node-date intervals. <a href="{displayed_timetree}">Open the full-resolution dated phylogeny</a>.</figcaption></figure>
            <details class="evidence-files" open><summary>Evidence and downloads</summary>{confidence_downloads}</details></div>
          </div>
        </section>
        """
    else:
        timetree_visual = f"""
        <section class="stage" id="time-tree" data-stage="3">
          <div class="stage-index"><span>3</span><b>Time scale</b></div>
          <div class="stage-body">
            <div class="stage-head">
              <div><h2>Time scaling was not performed</h2><p class="question">The workflow stops here when the temporal-signal gate does not pass.</p></div>
              <strong class="decision stop">DO NOT TIME-SCALE</strong>
            </div>
            <p>{escape(str(assessment["reason"]))}</p>
            <p>An undated recombination-corrected phylogeny, topology, and pairwise clonal SNP
            evidence remain available. Calendar dates for internal nodes would not be defensible
            from this analysis.</p>
          </div>
        </section>
        """

    outliers = directory / "clock" / "outliers.tsv"
    outlier_note = (
        "Automatic clock-outlier filtering was disabled for the temporal-signal gate, so all "
        "supplied genomes contributed. Review clock deviations in the root-to-tip table."
    )
    if outliers.is_file():
        rows = [line for line in outliers.read_text(encoding="utf-8").splitlines() if line]
        count = max(0, len(rows) - 1)
        outlier_note = (
            f"TreeTime reported {count} potential clock outlier(s). Review the outlier table."
        )

    raw_nearest = context.get("nearest_screening_contexts", [])
    nearest = raw_nearest if isinstance(raw_nearest, list) else []
    nearest_rows: list[str] = []
    for value in nearest:
        if not isinstance(value, dict):
            continue
        distance = value.get("min_ska_distance")
        distance_text = "—" if distance is None else f"{float(distance):.3g}"
        nearest_rows.append(
            "<tr>"
            f"<td>{escape(str(value.get('sample_id', '')))}</td>"
            f"<td>{escape(str(value.get('country', '') or 'Unknown'))}</td>"
            f"<td>{escape(str(value.get('collection_date', '') or 'Unknown'))}</td>"
            f"<td>{escape(str(value.get('nearest_focal', '')))}</td>"
            f"<td>{escape(distance_text)}</td>"
            "</tr>"
        )
    if context_count:
        nearest_table = (
            _table_region(
                "<table><thead><tr><th>Context genome</th><th>Country</th><th>Date</th>"
                "<th>Nearest focal isolate in screen</th><th>SKA screening distance</th>"
                "</tr></thead><tbody>" + "".join(nearest_rows) + "</tbody></table>",
                label="Public contextual genomes",
            )
            if nearest_rows
            else "<p>The context set was supplied without a ChronoClade selection manifest.</p>"
        )
        context_visual = f"""
        <section class="card">
          <h3>Public contextual genomes</h3>
          <p>{escape(str(context.get("interpretation", "")))}</p>
          <p><strong>{context_count}</strong> context genomes cover <strong>{len(context_locations)}</strong>
          reported locations. These are the nearest genomes within the downloaded screening pool,
          not necessarily the nearest in all public data. SKA distances were used only to select
          candidates; final relatedness must be interpreted from the recombination-corrected tree.</p>
          {nearest_table}
        </section>
        """
    else:
        context_visual = f"""
        <section class="card caveat">
          <h3>Public contextual genomes</h3>
          <p>{escape(str(context.get("interpretation", "No contextual genomes were supplied.")))}</p>
          <p>The report can describe structure among local isolates, but it cannot distinguish
          persistence from repeated introductions without an appropriate same-lineage context set.</p>
        </section>
        """

    randomised_values = [
        value for value in temporal.get("randomised", []) if isinstance(value, dict)
    ]
    exceedances = sum(
        float(value["r_squared"]) >= observed_r_squared for value in randomised_values
    )
    root_decision = "CONTINUE TO TEST" if observed_rate > 0 else "STOP"
    root_class = "review" if observed_rate > 0 else "stop"
    final_decision = "PROCEED" if bool(assessment["supported"]) else "DO NOT TIME-SCALE"
    final_class = "proceed" if bool(assessment["supported"]) else "stop"
    rtt_downloads = _download_list(
        [
            (
                "clock/root_to_tip_regression.svg",
                "Root-to-tip figure",
                "TreeTime scalable vector figure",
            ),
            ("root_to_tip.png", "Root-to-tip figure", "High-resolution raster figure"),
            ("clock/rtt.csv", "Root-to-tip data", "Collection dates, divergence and residuals"),
            ("clock/molecular_clock.txt", "Exploratory clock fit", "Original TreeTime statistics"),
        ],
        directory,
    )
    randomisation_downloads = _download_list(
        [
            ("date_randomisation.svg", "Randomisation figure", "Scalable vector figure"),
            ("date_randomisation.png", "Randomisation figure", "High-resolution raster figure"),
            (
                "date_randomisation.csv",
                "All randomisation results",
                "Observed and every permuted fit",
            ),
            (
                "temporal_signal.json",
                "Temporal-signal record",
                "Machine-readable gate input and result",
            ),
        ],
        directory,
    )
    audit_downloads = _download_list(
        [
            (
                "clonalframeml.labelled_tree.newick",
                "Corrected phylogeny",
                "Recombination-corrected Newick tree",
            ),
            (
                "clonalframeml.importation_status.txt",
                "Recombination calls",
                "ClonalFrameML importation output",
            ),
            ("clonal_snp_heatmap.svg", "Distance heatmap", "Scalable vector figure"),
            ("clonal_snp_heatmap.png", "Distance heatmap", "High-resolution raster figure"),
            (
                "clonal_pairwise_distances.tsv",
                "Pairwise distance table",
                "Clonal SNPs and callable sites",
            ),
            ("clonal_snp_matrix.tsv", "Clonal SNP matrix", "Full square matrix"),
            ("pairwise_callable_sites.tsv", "Callable-site matrix", "Denominators for every pair"),
            (
                "public_health_evidence.json",
                "Evidence ledger",
                "Machine-readable interpretation inputs",
            ),
            ("report.json", "Analysis record", "Machine-readable lineage report"),
            ("context_manifest.tsv", "Context manifest", "Frozen public-context selection"),
        ],
        directory,
    )
    write_supporting_bundle(directory)
    font_css = _embedded_font_css()
    species_display = str(report["species"]).replace("_", " ")
    lineage_display = str(report["lineage"])
    output = directory / "report.html"
    output.write_text(
        f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ChronoClade analysis report — {escape(species_display)} {escape(lineage_display)}</title>
  <style>
    {font_css}
    :root {{ color-scheme:light; --ink:#17191f; --muted:#5f6470; --line:#c9ccd4; --paper:#fff; --wash:#eef0f4; --blue:#2855a6; --red:#d63c2f; --amber:#9a6500; --green:#14734f; --rail:248px; }}
    * {{ box-sizing:border-box; }}
    html {{ scroll-behavior:smooth; }}
    body {{ margin:0; background:var(--wash); color:var(--ink); font:16px/1.58 Archivo, "Helvetica Neue", sans-serif; }}
    ::selection {{ background:#ffdf78; color:var(--ink); }}
    a {{ color:var(--blue); text-underline-offset:3px; }} a:focus-visible,summary:focus-visible {{ outline:3px solid #ffb800; outline-offset:3px; }}
    .shell {{ width:min(1480px,100%); margin:auto; background:var(--paper); min-height:100vh; }}
    .identity {{ display:grid; grid-template-columns:var(--rail) minmax(0,1fr); border-bottom:2px solid var(--ink); }}
    .identity-mark {{ background:var(--ink); color:#fff; padding:30px 24px; font-weight:800; letter-spacing:.06em; }}
    .identity-copy {{ padding:24px 42px 26px; display:flex; flex-direction:column; gap:14px; align-items:flex-start; min-width:0; }}
    h1 {{ margin:0; font-size:clamp(32px,5vw,68px); line-height:1.12; letter-spacing:-.035em; max-width:900px; }} h1 i {{ font-style:italic; }} h1 span {{ display:block; }}
    .identity-copy p {{ margin:10px 0 0; color:var(--muted); max-width:70ch; }}
    .overall {{ width:max-content; max-width:100%; border:1px solid var(--red); padding:9px 12px; color:var(--red); }} .overall.supported {{ border-color:var(--red); color:var(--red); }} .overall small {{ color:currentColor; }}
    .overall small,.measure-strip small {{ display:block; text-transform:uppercase; letter-spacing:.08em; color:var(--muted); font-size:11px; font-weight:700; }}
    .overall b {{ display:block; font-size:21px; margin-top:4px; }}
    .contents {{ position:sticky; top:0; z-index:5; display:flex; gap:0; padding-left:var(--rail); background:#fff; border-bottom:1px solid var(--ink); overflow:auto; }}
    .contents a {{ flex:1; min-width:170px; padding:13px 18px; border-left:1px solid var(--line); color:var(--ink); text-decoration:none; font-size:13px; font-weight:700; }} .contents a:hover {{ background:#edf3ff; }}
    .stage {{ display:grid; grid-template-columns:var(--rail) minmax(0,1fr); border-bottom:2px solid var(--ink); scroll-margin-top:52px; }}
    .stage-index {{ padding:38px 24px; border-right:1px solid var(--ink); background:#f4f5f7; }}
    .stage-index span {{ display:block; font-size:74px; line-height:.8; font-weight:800; letter-spacing:-.06em; color:var(--blue); }}
    .stage-index b {{ display:block; margin-top:20px; text-transform:uppercase; letter-spacing:.1em; font-size:12px; }}
    .stage.active .stage-index {{ background:var(--ink); color:#fff; }} .stage.active .stage-index span {{ color:#fff; }}
    .stage-body {{ padding:40px 44px 50px; min-width:0; }}
    .stage-head {{ display:flex; justify-content:space-between; gap:28px; align-items:start; margin-bottom:22px; }} .stage-head>div {{ min-width:0; }}
    h2 {{ margin:0; max-width:24ch; font-size:clamp(27px,3.2vw,44px); line-height:1.05; letter-spacing:-.025em; overflow-wrap:break-word; text-wrap:balance; }}
    h3 {{ margin:42px 0 8px; font-size:21px; }} p {{ max-width:75ch; }} .question {{ color:var(--muted); font-size:18px; margin:8px 0 0; }}
    .decision {{ display:block; min-width:150px; padding:10px 12px; border:2px solid currentColor; text-align:center; font-size:12px; letter-spacing:.08em; }}
    .decision.proceed,.decision.stop {{ color:var(--red); }} .decision.review {{ color:var(--blue); }}
    figure {{ margin:28px 0 0; min-width:0; }} figure img {{ display:block; width:100%; height:auto; border:1px solid var(--ink); background:#fff; }} figcaption {{ margin-top:8px; color:var(--muted); font-size:13px; }}
    .evidence-layout {{ display:grid; grid-template-columns:minmax(0,1fr) 270px; gap:22px; align-items:start; min-width:0; }} .evidence-layout>* {{ min-width:0; }} .evidence-layout .evidence-files {{ grid-column:2; grid-row:1/4; margin-top:28px; min-width:0; }}
    .evidence-layout .download a {{ grid-template-columns:minmax(0,1fr) 42px; }} .evidence-layout .download small {{ grid-column:1/-1; grid-row:2; }}
    .figure-scroll {{ max-width:100%; overflow:auto; scrollbar-color:var(--blue) #e5e7ec; }} .figure-scroll:focus-visible,.table-scroll:focus-visible {{ outline:3px solid #ffb800; outline-offset:3px; }}
    .measure-strip {{ display:grid; grid-template-columns:repeat(3,1fr); border:1px solid var(--ink); margin:24px 0; }}
    .measure-strip > div {{ padding:16px; border-right:1px solid var(--line); }} .measure-strip > div:last-child {{ border:0; }}
    .measure-strip b {{ display:block; font-size:21px; margin:4px 0 2px; font-variant-numeric:tabular-nums; }} .measure-strip span {{ color:var(--muted); font-size:12px; }}
    .confidence-grid {{ display:grid; grid-template-columns:repeat(3,1fr); border:1px solid var(--ink); margin:24px 0; }} .confidence-grid>div {{ padding:16px; border-right:1px solid var(--line); border-bottom:1px solid var(--line); }} .confidence-grid>div:nth-child(3n) {{ border-right:0; }} .confidence-grid>div:nth-last-child(-n+3) {{ border-bottom:0; }} .confidence-grid small {{ display:block; text-transform:uppercase; letter-spacing:.08em; color:var(--muted); font-size:11px; font-weight:700; }} .confidence-grid b {{ display:block; font-size:21px; margin:4px 0 2px; font-variant-numeric:tabular-nums; }} .confidence-grid span {{ color:var(--muted); font-size:12px; }}
    .logic {{ display:grid; grid-template-columns:1fr auto 1fr; align-items:center; gap:18px; margin:26px 0; }} .logic div {{ padding:20px; background:#f4f5f7; border:1px solid var(--line); }} .logic b {{ display:block; font-size:18px; }} .logic span {{ font-size:30px; color:var(--muted); }}
    details.evidence-files {{ margin-top:22px; border-top:1px solid var(--ink); border-bottom:1px solid var(--ink); }} details.evidence-files summary {{ padding:14px 0; cursor:pointer; font-weight:800; }}
    .downloads {{ list-style:none; margin:0 0 14px; padding:8px 12px 12px; border-top:1px solid var(--line); }} .download a {{ display:grid; grid-template-columns:minmax(180px,1fr) 2fr 54px; gap:14px; padding:11px 4px; border-bottom:1px solid var(--line); text-decoration:none; align-items:center; }} .download small {{ color:var(--muted); }} .download b {{ text-align:right; font-size:11px; letter-spacing:.08em; color:var(--muted); }}
    .scenario,.card,.verdict {{ margin:22px 0 0; padding:22px 0; border-top:1px solid var(--ink); background:#fff; }} .scenario strong {{ display:block; font-size:28px; max-width:34ch; }} .confidence,.evidence {{ font-weight:800; }} .guardrail,.caveat {{ background:#fff8dc; padding:15px; }} .card img {{ display:block; max-width:100%; height:auto; }}
    table {{ width:100%; border-collapse:collapse; margin-top:16px; font-size:14px; }} th,td {{ padding:10px 8px; text-align:left; border-bottom:1px solid var(--line); vertical-align:top; }} th {{ font-size:11px; text-transform:uppercase; letter-spacing:.06em; color:var(--muted); }}
    .table-scroll {{ max-width:100%; overflow:auto; scrollbar-color:var(--blue) #e5e7ec; }} .scroll-hint {{ display:none; }}
    .report-close {{ display:grid; grid-template-columns:var(--rail) 1fr; background:var(--ink); color:#fff; }} .report-close b {{ padding:30px 24px; border-right:1px solid #555b67; }} .report-close div {{ padding:30px 42px; }} .report-close a {{ color:#fff; }}
    #root-to-tip .stage-body {{ display:grid; grid-template-columns:minmax(270px,.72fr) minmax(560px,1.55fr); column-gap:30px; align-items:start; }} #root-to-tip .stage-head,#root-to-tip .stage-body>p,#root-to-tip .measure-strip {{ grid-column:1; }} #root-to-tip .stage-head {{ display:block; }} #root-to-tip .decision {{ margin-top:18px; width:max-content; }} #root-to-tip .measure-strip {{ grid-template-columns:1fr; }} #root-to-tip .measure-strip>div {{ border-right:0; border-bottom:1px solid var(--line); }} #root-to-tip .evidence-layout {{ grid-column:2; grid-row:1/5; grid-template-columns:minmax(0,1fr) 210px; }} #root-to-tip .logic {{ grid-column:1; }}
    @media(max-width:1100px) {{ #root-to-tip .stage-body {{ display:block; }} .evidence-layout,#root-to-tip .evidence-layout {{ grid-template-columns:1fr; }} .evidence-layout .evidence-files {{ grid-column:1; grid-row:auto; margin-top:8px; }} }}
    @media(max-width:800px) {{ :root {{ --rail:76px; }} .identity-copy {{ display:block; padding:22px 20px; }} .overall {{ margin-top:22px; }} .contents {{ padding-left:0; overflow:visible; }} .contents a {{ min-width:0; padding:12px 7px; text-align:center; font-size:11px; white-space:nowrap; }} .stage-index {{ padding:28px 10px; }} .stage-index span {{ font-size:48px; }} .stage-index b {{ writing-mode:vertical-rl; margin:18px auto 0; }} .stage-body {{ padding:30px 18px 38px; }} .stage-head {{ display:block; }} .decision {{ margin-top:18px; width:max-content; max-width:100%; }} .measure-strip,.confidence-grid {{ grid-template-columns:1fr; }} .measure-strip > div,.confidence-grid>div {{ border-right:0; border-bottom:1px solid var(--line); min-width:0; }} .confidence-grid b {{ overflow-wrap:anywhere; }} .confidence-grid>div:nth-child(3n) {{ border-right:0; }} .confidence-grid>div:nth-last-child(-n+3) {{ border-bottom:1px solid var(--line); }} .confidence-grid>div:last-child {{ border-bottom:0; }} .logic {{ grid-template-columns:1fr; }} .logic > span {{ transform:rotate(90deg); justify-self:center; }} .download a {{ grid-template-columns:minmax(0,1fr) 42px; }} .download span,.download small {{ overflow-wrap:anywhere; }} .download small {{ grid-column:1/-1; grid-row:2; }} .scroll-hint {{ display:block; position:sticky; left:0; width:max-content; max-width:100%; margin-top:12px; padding:9px 10px; background:#edf3ff; color:var(--blue); font-size:12px; font-weight:700; }} .table-scroll:not(.compact-table) table {{ min-width:680px; }} .table-scroll:not(.compact-table) th:first-child,.table-scroll:not(.compact-table) td:first-child {{ position:sticky; left:0; background:#fff; z-index:1; }} .compact-table .scroll-hint,.compact-table thead {{ display:none; }} .compact-table table,.compact-table tbody,.compact-table tr,.compact-table td {{ display:block; min-width:0; width:100%; }} .compact-table tr {{ padding:10px 0; border-bottom:1px solid var(--line); }} .compact-table td {{ padding:3px 0; border:0; }} .compact-table td:first-child {{ font-weight:700; }} .figure-scroll img {{ min-width:760px; }} }}
    @media(prefers-reduced-motion:reduce) {{ html {{ scroll-behavior:auto; }} }}
    @media print {{ .contents {{ display:none; }} body,.shell {{ background:#fff; }} .stage {{ break-inside:avoid-page; }} details {{ display:block; }} details > * {{ display:block; }} .report-close {{ color:#000; background:#fff; border-top:2px solid #000; }} .report-close a {{ color:#000; }} }}
  </style>
</head>
<body>
<main class="shell">
  <header class="identity">
    <div class="identity-mark">ChronoClade<br>ANALYSIS</div>
    <div class="identity-copy"><div><h1><i>{escape(species_display)}</i><span>{escape(lineage_display)}</span></h1><p>{int(report["sample_count"])} genomes · {int(report["distinct_dates"])} distinct collection dates · recombination-aware temporal analysis</p></div><div class="overall {escape(status)}"><small>Temporal decision</small><b>{escape(str(assessment["label"]))}</b></div></div>
  </header>
  <nav class="contents" aria-label="Analysis stages"><a href="#root-to-tip">1 · Explore</a><a href="#randomisation">2 · Test</a><a href="#time-tree">3 · Time scale</a><a href="#interpretation">4 · Interpret</a></nav>

  <section class="stage active" id="root-to-tip" data-stage="1">
    <div class="stage-index"><span>1</span><b>Explore</b></div>
    <div class="stage-body">
      <div class="stage-head"><div><h2>Does divergence increase with sampling time?</h2><p class="question">Root-to-tip regression is the first diagnostic—not the final proof.</p></div><strong class="decision {root_class}">{root_decision}</strong></div>
      <p>Each point is one genome. Its vertical position is genetic distance from the fitted root; its horizontal position is collection date. A positive slope and a coherent cloud suggest that the dates may contain clock information. The plot can also reveal extreme residuals or obvious structure.</p>
      <div class="measure-strip"><div><small>Exploratory rate</small><b>{observed_rate:.3g}</b><span>substitutions/site/year</span></div><div><small>Exploratory root-to-tip fit</small><b>R² {observed_r_squared:.3g}</b><span>descriptive, not a pass criterion alone</span></div><div><small>Samples</small><b>{int(report["sample_count"])}</b><span>{int(report["distinct_dates"])} distinct dates</span></div></div>
      <div class="evidence-layout"><figure>{root_visual}<figcaption>Root-to-tip regression after recombination correction and least-squares rooting. {escape(outlier_note)}</figcaption></figure>
      <div class="logic"><div><b>What this can show</b>Genetic divergence is associated with collection time.</div><span>→</span><div><b>What it cannot show alone</b>That the association is stronger than a result obtained with arbitrary dates.</div></div>
      <details class="evidence-files" open><summary>Evidence and downloads</summary>{rtt_downloads}</details></div>
    </div>
  </section>

  <section class="stage" id="randomisation" data-stage="2">
    <div class="stage-index"><span>2</span><b>Test</b></div>
    <div class="stage-body">
      <div class="stage-head"><div><h2>Is the observed fit stronger than shuffled dates?</h2><p class="question">This is the formal proceed/stop gate for time scaling.</p></div><strong class="decision {final_class}">{final_decision}</strong></div>
      <p>The same corrected tree is analysed repeatedly after collection dates are permuted among genomes. The empirical p-value is the proportion of shuffled analyses whose root-to-tip R² equals or exceeds the observed value, with a one-count correction.</p>
      <div class="measure-strip"><div><small>Exploratory root-to-tip fit</small><b>R² {observed_r_squared:.3g}</b><span>real collection dates</span></div><div><small>Null exceedances</small><b>{exceedances} / {successful}</b><span>successful randomisations</span></div><div><small>Empirical p</small><b>{p_value:.3g}</b><span>predefined threshold {p_value_threshold:.3g}</span></div></div>
      <div class="evidence-layout"><figure><img src="date_randomisation.svg" alt="Histograms comparing observed TreeTime fit and rate with date-randomised values"><figcaption>The red line is the observed result; blue bars are the null distribution generated by permuting collection dates.</figcaption></figure>
      <section class="verdict {escape(status)}"><strong>{escape(str(assessment["label"]))}</strong><p>{escape(str(assessment["reason"]))}</p></section>
      <details class="evidence-files" open><summary>Evidence and downloads</summary>{randomisation_downloads}</details></div>
    </div>
  </section>

  {timetree_visual}

  <section class="stage" id="interpretation" data-stage="4">
    <div class="stage-index"><span>4</span><b>Interpret</b></div>
    <div class="stage-body">
      <div class="stage-head"><div><h2>What pattern is consistent with these genomes?</h2><p class="question">Time, topology, distance and context are considered together.</p></div></div>
      {public_health_summary}
      {context_visual}
      {public_health_evidence}
      <section class="card caveat"><h3>Interpretation boundary</h3><p>Temporal signal means the dates inform evolutionary rate and node timing under this model. It does not prove direct transmission, local circulation, or a definitive number of introductions. Those conclusions also depend on recombination filtering, contextual sampling and epidemiological metadata.</p></section>
      <h3>Complete audit package</h3><p>Every reader-facing result is available separately and as one ZIP bundle.</p>{audit_downloads}
      <p><a href="supporting_results.zip"><strong>Download all supporting results (.zip)</strong></a></p>
    </div>
  </section>

  <footer class="report-close"><b>REPORT COMPLETE</b><div><p>Generated by ChronoClade. The dated-tree figure, node-date intervals, clock statistics and source tables are preserved in the supporting-results bundle.</p></div></footer>
</main>
<script>
  const stages = [...document.querySelectorAll('.stage')];
  const observer = new IntersectionObserver((entries) => {{
    const visible = entries.filter(entry => entry.isIntersecting)
      .sort((a, b) => b.intersectionRatio - a.intersectionRatio)[0];
    if (!visible) return;
    stages.forEach(stage => stage.classList.toggle('active', stage === visible.target));
  }}, {{ rootMargin: '-18% 0px -58% 0px', threshold: [0, .1, .35, .6] }});
  stages.forEach(stage => observer.observe(stage));
</script>
</body>
</html>
""",
        encoding="utf-8",
    )
    return output


def write_summary_report(summary: dict[str, object], *, output: Path) -> Path:
    """Write a compact landing page linking all lineage reports."""

    rows: list[str] = []
    lineages = summary.get("lineages", [])
    if not isinstance(lineages, list):
        raise ValueError("summary lineages must be a list")
    for lineage in lineages:
        if not isinstance(lineage, dict):
            continue
        status = str(lineage.get("temporal_status", lineage.get("status", "unknown")))
        label = status.replace("_", " ").title()
        metrics = lineage.get("temporal_signal", {})
        rate = "—"
        r_squared = "—"
        p_value = "—"
        if isinstance(metrics, dict) and "observed" in metrics:
            rate = _format_rate(_metric(metrics, "rate"))
            r_squared = f"{_metric(metrics, 'r_squared'):.3g}"
            p_value = f"{float(metrics['p_value_r_squared']):.3g}"
        slug = str(lineage["slug"])
        lineage_context = lineage.get("context", {})
        context_count = (
            int(lineage_context.get("context_samples", 0))
            if isinstance(lineage_context, dict)
            else 0
        )
        lineage_public_health = lineage.get("public_health", {})
        public_health = lineage_public_health if isinstance(lineage_public_health, dict) else {}
        lineage_scenario = public_health.get("scenario", {})
        scenario = lineage_scenario if isinstance(lineage_scenario, dict) else {}
        scenario_label = str(scenario.get("label", "Not assessed"))
        scenario_code = str(scenario.get("code", "indeterminate"))
        report_link = (
            f'<a href="{escape(slug, quote=True)}/report.html">Open report</a>'
            if "temporal_signal" in lineage
            else "Not analysed"
        )
        rows.append(
            f"<tr><td>{escape(str(lineage['species']))}</td><td>{escape(str(lineage['lineage']))}</td>"
            f"<td>{int(lineage['sample_count'])}</td><td>{context_count}</td>"
            f'<td><span class="status {escape(scenario_code)}">{escape(scenario_label)}</span></td>'
            f'<td><span class="status {escape(status)}">{escape(label)}</span></td>'
            f"<td>{escape(rate)}</td><td>{escape(r_squared)}</td><td>{escape(p_value)}</td><td>{report_link}</td></tr>"
        )

    report_path = output / "index.html"
    report_path.write_text(
        f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>ChronoClade analysis report</title>
<style>
body{{margin:0;background:#f2f6f8;color:#1d2b34;font:16px/1.5 Inter,ui-sans-serif,system-ui,sans-serif}}main{{width:min(1200px,calc(100% - 32px));margin:32px auto}}header{{background:#17394d;color:white;padding:30px;border-radius:16px}}h1{{margin:0}}header p{{color:#d5e5ed}}section{{background:white;border:1px solid #d8e1e7;border-radius:13px;margin-top:18px;padding:22px;overflow:auto}}table{{width:100%;border-collapse:collapse;min-width:900px}}th,td{{padding:12px;text-align:left;border-bottom:1px solid #e4eaee}}th{{font-size:13px;text-transform:uppercase;color:#5b6b76}}a{{color:#265f7d;font-weight:700}}.status{{font-weight:700}}.supported{{color:#157347}}.not_supported{{color:#a33a36}}.not_assessed{{color:#8a6420}}footer{{margin-top:20px;color:#5b6b76}}
</style></head><body><main><header><h1>ChronoClade analysis report</h1><p>Recombination-aware temporal and contextual bacterial phylogenetics</p></header>
<section><h2>Lineage results</h2><table><thead><tr><th>Species</th><th>Lineage</th><th>Genomes</th><th>Context</th><th>Working interpretation</th><th>Temporal result</th><th>Rate</th><th>R²</th><th>p</th><th></th></tr></thead><tbody>{"".join(rows)}</tbody></table></section>
<section><h2>How to interpret this report</h2><p>The working interpretation combines recombination-filtered distances, corrected topology, longitudinal sampling and public context. It is provisional and must be reviewed with epidemiological information. A supported temporal signal permits time scaling under the fitted model but is not required for the genomic scenario assessment.</p></section>
<footer>Generated by ChronoClade.</footer></main></body></html>
""",
        encoding="utf-8",
    )
    return report_path
