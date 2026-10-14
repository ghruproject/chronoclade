"""Presentation of saved LIN-group calibration results; no clock fitting."""

import csv
from html import escape

import numpy as np

from chronoclade.report_components.stage import write_report

NOTE = (
    "LIN6 and LIN7 groups are nested analyses, not independent confirmations. "
    "All genomes in each matching frozen group contribute genetics to its NJ tree. "
    "Collection dates from 1900 onward are eligible for regression only; metadata and "
    "genetic membership are retained. Each latest-year fold chooses its root and line "
    "using strictly earlier collection dates, then predicts root-to-tip allele distances from withheld dates. "
    "This is a transductive held-out-date check, not validation on unseen genomes. "
    "Root choice uses dates, tips share ancestry, and no temporal significance or "
    "independent-tip confidence intervals are claimed."
)


def _plots(group, path):
    from chronoclade.plotting import headless_pyplot

    plt = headless_pyplot()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    fit = group["full_date_fit"]
    points = fit["points"]
    x = np.asarray([p["midpoint_year"] for p in points])
    axes[0].scatter(x, [p["root_to_tip"] for p in points], s=6, alpha=0.25, color="#686174")
    focal = [p for p in points if p.get("role") == "input"]
    axes[0].scatter(
        [p["midpoint_year"] for p in focal],
        [p["root_to_tip"] for p in focal],
        s=14,
        color="#773886",
        label="Input genomes",
    )
    ends = np.array([min(x), max(x)])
    axes[0].plot(
        ends,
        fit["centered_intercept"] + fit["slope"] * (ends - fit["reference_year"]),
        color="#242026",
    )
    axes[0].set(
        xlabel="Collection-date midpoint year",
        ylabel="Estimated allele-distance units",
        title=f"Full group · rate {fit['slope']:.3g}/year · R² {fit['r_squared'] or 0:.3g}",
    )
    axes[0].legend(frameon=False)
    folds = group["heldout_dates"]
    evaluated = [f for f in folds if f["model_distance_mae"] is not None]
    if evaluated:
        positions = np.arange(len(evaluated))
        axes[1].bar(
            positions - 0.18,
            [f["model_distance_mae"] for f in evaluated],
            0.36,
            color="#773886",
            label="Training-root accumulation prediction",
        )
        axes[1].bar(
            positions + 0.18,
            [f["baseline_distance_mae"] for f in evaluated],
            0.36,
            color="#bbb4c3",
            label="Same-root training-mean distance",
        )
        axes[1].set_xticks(positions, [f["held_out_year"] for f in evaluated])
        axes[1].legend(frameon=False, fontsize=8)
    else:
        axes[1].text(
            0.5,
            0.5,
            "No evaluable forward-prediction folds",
            ha="center",
            transform=axes[1].transAxes,
        )
    axes[1].set(
        xlabel="Held-out calendar-year block",
        ylabel="Mean absolute distance error (allele units)",
        title="Transductive held-out-date check",
    )
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def _group_label(group):
    return f"LIN{group['level']} " + ".".join(map(str, group["prefix"]))


def _evaluated(group):
    return [
        fold
        for fold in group["heldout_dates"]
        if all(
            fold.get(field) is not None
            for field in ("model_distance_mae", "baseline_distance_mae", "persistence_distance_mae")
        )
    ]


def _beats_both(fold):
    return fold["model_distance_mae"] < min(
        fold["baseline_distance_mae"], fold["persistence_distance_mae"]
    )


def _number(value, precision=3):
    return f"{value:.{precision}g}" if value is not None else "unavailable"


def _fold_details(group):
    rows = []
    for fold in group["heldout_dates"]:
        cells = [
            fold["held_out_year"],
            fold["n_training"],
            fold["n_test"],
            _number(fold.get("training_rate")),
            _number(fold["model_distance_mae"]),
            _number(fold["baseline_distance_mae"]),
            _number(fold["persistence_distance_mae"]),
            fold["reason"]
            or ("Beat both baselines" if _beats_both(fold) else "Did not beat both baselines"),
        ]
        rows.append(
            "<tr>" + "".join("<td>" + escape(str(cell)) + "</td>" for cell in cells) + "</tr>"
        )
    if not rows:
        return "<p>" + escape(group["reason"] or "No held-out-year checks available.") + "</p>"
    return (
        '<div class="table-scroll"><table><thead><tr><th>Held-out year</th>'
        "<th>Training genomes</th><th>Test genomes</th><th>Training rate / year</th>"
        "<th>Model MAE</th><th>Mean baseline MAE</th><th>Persistence MAE</th>"
        "<th>Descriptive result</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table></div>"
        '<p class="muted">All errors are estimated allele-distance units under the same frozen '
        "training root. Collection-date inversion remains a secondary saved diagnostic.</p>"
    )


def _report(result, output):
    """Render saved estimates and errors without recomputing trees or fits."""
    flat, figures = [], []
    folds = [fold for group in result["groups"] for fold in _evaluated(group)]
    wins = sum(_beats_both(fold) for fold in folds)
    headline = f"Model beat both baselines in {wins} of {len(folds)} held-out group/year checks."
    fragments = [
        "<p>Rates are estimated allele-distance units/year. These nested groups share "
        "genomes; their checks are not independent confirmations.</p>",
        '<details><summary>How the prediction check works</summary>',
        "<p>" + escape(NOTE) + "</p>",
    ]
    loci = sorted(
        {
            group["normalization_loci"]
            for group in result["groups"]
            if group.get("normalization_loci")
        }
    )
    fragments.append(
        "<p>Rates are estimated allele-distance units/year, scaled by recorded locus "
        "counts ("
        + escape(", ".join(map(str, loci)))
        + " loci in the frozen export), not literal allele "
        "events. Pair-specific callable denominators vary. Forward distance prediction uses "
        "training dates only and is compared with same-root training-mean distance and "
        "latest-training-year mean-distance persistence baselines. No biological cutoff is recommended.</p></details>"
    )
    fragments.append(
        '<div class="table-scroll"><table><thead><tr><th>Group</th>'
        "<th>Inputs / public</th><th>Dated / all</th><th>Rate / year</th><th>R²</th>"
        "<th>Model beats both baselines</th></tr></thead><tbody>"
    )
    for group in result["groups"]:
        fit = group["full_date_fit"] or {}
        checks = group["heldout_dates"]
        warnings = [
            f"{fold['held_out_year']}: {fold['reason']}" for fold in checks if fold["reason"]
        ]
        warnings += [
            f"{fold['held_out_year']}: distance MAE exceeds mean baseline"
            for fold in checks
            if fold["model_distance_mae"] is not None
            and fold["model_distance_mae"] > fold["baseline_distance_mae"]
        ]
        summary = "; ".join(
            f"{fold['held_out_year']}: {fold['n_training']}/{fold['n_test']}; "
            + " / ".join(
                _number(fold[field], 5)
                for field in (
                    "model_distance_mae",
                    "baseline_distance_mae",
                    "persistence_distance_mae",
                )
            )
            for fold in checks
        )
        evaluated = _evaluated(group)
        group_wins = sum(_beats_both(fold) for fold in evaluated)
        row = {
            "group": _group_label(group),
            "inputs": group["n_inputs"],
            "public": group["n_public"],
            "dated": group["n_dated"],
            "all_genomes": group["n_genomes"],
            "rate": fit.get("slope"),
            "r_squared": fit.get("r_squared"),
            "status": group["status"],
            "heldout_summary": summary,
            "warnings": "; ".join(warnings)
            or group["reason"]
            or "Descriptive; see fold sensitivity",
            "wins_against_both_baselines": group_wins,
            "evaluated_heldout_years": len(evaluated),
        }
        flat.append(row)
        cells = [
            row["group"],
            f"{row['inputs']} / {row['public']}",
            f"{row['dated']} / {row['all_genomes']}",
            _number(row["rate"]),
            _number(row["r_squared"]),
            f"{group_wins} / {len(evaluated)} years" if evaluated else "insufficient dates",
        ]
        fragments.append(
            "<tr>" + "".join("<td>" + escape(str(cell)) + "</td>" for cell in cells) + "</tr>"
        )
        if group["full_date_fit"]:
            figures.append(
                (
                    group["group_id"] + "/calibration.svg",
                    row["group"] + " · full group and held-out dates",
                )
            )
    fragments.append("</tbody></table></div>")
    for group in result["groups"]:
        fragments.append(
            "<details><summary>"
            + escape(_group_label(group) + " · year checks, inputs and uncertainty")
            + "</summary>"
            + _fold_details(group)
            + "<p>Inputs: "
            + escape(", ".join(group["input_accessions"]))
            + "</p><p>"
            + escape(group.get("uncertainty_note", group["reason"] or ""))
            + "</p><p>"
            + escape(
                "Training-fold rate range: "
                + str(group.get("training_fold_rate_range"))
                + "; not a confidence interval. Distinct training-fold root edges: "
                + str(group.get("training_fold_distinct_root_edge_count"))
            )
            + "</p></details>"
        )
    if flat:
        with (output / "summary.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(flat[0]))
            writer.writeheader()
            writer.writerows(flat)
    report = dict(
        result,
        selected_sample_ids=[],
        temporal_assessment={
            "code": "sensitivity",
            "supported": False,
            "reason": headline
            + " Nested group/year checks are not independent; empirical distance cutoffs are not established.",
        },
    )
    return write_report(
        output / "report.html",
        title="LIN6 / LIN7 temporal calibration trial",
        result=report,
        overview=f"{len(flat)} input-containing frozen groups · all group genomes retained",
        report_label="TEMPORAL CALIBRATION TRIAL",
        figures=figures,
        extra="".join(fragments),
        figure_heading="Accumulation and held-out-year checks",
        assessment_first=True,
        assessment_title="Temporal calibration assessment",
    )
