"""Optional independent branch-support inference on the filtered alignment."""

import csv
from html import escape
from pathlib import Path

from Bio import Phylo

from chronoclade.errors import WorkflowError
from chronoclade.stage_artifacts import reference


def validate_replicates(replicates):
    if isinstance(replicates, bool) or not isinstance(replicates, int) or (replicates != 0 and replicates < 1000):
        raise WorkflowError("Ultrafast bootstrap replicates must be 0 (disabled) or at least 1000")


def splits(tree):
    """Canonical non-trivial unrooted splits; root placement does not change keys."""
    tips = [tip.name for tip in tree.get_terminals()]
    if len(tips) != len(set(tips)) or any(not name for name in tips):
        raise WorkflowError("Branch-support tree has missing or duplicate tip identifiers")
    universe = frozenset(tips)
    result = {}
    for clade in tree.get_nonterminals():
        side = frozenset(tip.name for tip in clade.get_terminals())
        other = universe - side
        if min(len(side), len(other)) < 2:
            continue
        key = min(tuple(sorted(side)), tuple(sorted(other)), key=lambda x: (len(x), x))
        value = clade.confidence
        if value is not None and not 0 <= value <= 100:
            raise WorkflowError("Ultrafast bootstrap support is outside 0–100 percent")
        if key not in result or value is not None:
            result[key] = value
    return universe, result


def run_branch_support(*, alignment, corrected_tree, directory, output_root, replicates, threads, seed, force=False):
    """Infer a separate ML tree; never replace the ClonalFrameML dating topology."""
    from chronoclade.native_pipeline import _run_command

    validate_replicates(replicates)
    if not replicates:
        return {"status": "not_assessed", "reason": "Optional branch support disabled"}
    directory = Path(directory)
    prefix = directory / "filtered_support"
    supported_tree = prefix.with_suffix('.treefile')
    report = prefix.with_suffix('.iqtree')
    _run_command(
        ["iqtree", "-s", str(alignment), "-m", "GTR+G", "-B", str(replicates),
         "--bnni", "-seed", str(seed), "-T", str(threads), "-pre", str(prefix), "-redo"],
        log=directory / "logs" / "filtered_support.log", expected=(supported_tree, report),
        force=force, cwd=directory, inputs=(Path(alignment),),
    )
    original_ids, original = splits(Phylo.read(corrected_tree, 'newick'))
    support_ids, supported = splits(Phylo.read(supported_tree, 'newick'))
    if support_ids != original_ids:
        raise WorkflowError("Branch-support inference changed the exact selected tip set")
    comparison = directory / 'filtered_support_splits.csv'
    rows = [{"split_side": '|'.join(key), "present_in_support_tree": key in supported,
             "ufboot_percent": supported.get(key) if key in supported else ''} for key in sorted(original)]
    with comparison.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['split_side', 'present_in_support_tree', 'ufboot_percent'])
        writer.writeheader()
        writer.writerows(rows)
    return {
        "status": "complete", "method": "IQ-TREE UFBoot2 with bootstrap NNI", "replicates": replicates,
        "seed": seed, "model": "GTR+G", "alignment_basis": "ClonalFrameML-filtered alignment",
        "topology_policy": "Separate ML inference; corrected tree and dating topology unchanged",
        "corrected_splits": len(original), "shared_splits": len(set(original) & set(supported)),
        "shared_splits_with_support": sum(supported[key] is not None for key in set(original) & set(supported)),
        "inference_tree": reference(output_root, supported_tree),
        "iqtree_report": reference(output_root, report), "split_comparison": reference(output_root, comparison),
        "source_alignment": reference(output_root, alignment),
        "source_corrected_tree": reference(output_root, corrected_tree),
    }


def support_html(result):
    if result.get('status') != 'complete':
        return ''
    return (
        '<h3>Optional branch support</h3><p>'
        f"{escape(result['method'])}: {result['replicates']:,} replicates on the recombination-filtered alignment. "
        f"{result['shared_splits']} of {result['corrected_splits']} corrected-tree splits occur in the separate ML tree. "
        'Support percentages apply to matching splits only; an absent split has no assigned support. '
        'This inference does not replace the corrected tree or establish temporal signal. '
        'The supported Newick tree, IQ-TREE report and split comparison are in Downloads.</p>'
    )
