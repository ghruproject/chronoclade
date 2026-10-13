"""Project independently selected outgroup attachments onto an existing NJ tree."""
import copy


def project_outgroup_root(ingroup_tree, augmented_tree, outgroup_id):
    """Preserve ingroup distances/topology; transfer a supported basal split.

    The attachment fraction is estimated by NJ with the external outgroup, never
    fitted to dates. Incompatible topologies or negative basal lengths fail closed.
    """
    tree = copy.deepcopy(ingroup_tree)
    augmented = copy.deepcopy(augmented_tree)
    ingroup = {tip.name for tip in tree.get_terminals()}
    augmented_ids = {tip.name for tip in augmented.get_terminals()}
    if outgroup_id in ingroup or augmented_ids != ingroup | {outgroup_id}:
        raise ValueError('Root projection requires exactly the ingroup plus one external outgroup')
    augmented.root_with_outgroup(outgroup_id)
    tips = {tip.name: tip for tip in augmented.get_terminals()}
    ancestor = augmented.common_ancestor([tips[name] for name in ingroup])
    if {tip.name for tip in ancestor.get_terminals()} != ingroup or len(ancestor.clades) != 2:
        raise ValueError('Outgroup does not provide a resolved ingroup basal split')
    lengths = [float(child.branch_length or 0) for child in ancestor.clades]
    if min(lengths) < 0 or sum(lengths) <= 0:
        raise ValueError('Outgroup attachment has negative or unresolved basal branch lengths')
    sides = [{tip.name for tip in child.get_terminals()} for child in ancestor.clades]
    match = None
    for node in tree.find_clades():
        descendants = {tip.name for tip in node.get_terminals()}
        for number, side in enumerate(sides):
            if descendants == side:
                match = node, number
                break
        if match:
            break
    if match is None:
        raise ValueError('Outgroup-supported basal split is absent from the original ingroup tree')
    node, number = match
    original_length = float(node.branch_length or 0)
    if original_length <= 0:
        raise ValueError('Original ingroup root edge has nonpositive branch length')
    fraction = lengths[number] / sum(lengths)
    tree.root_with_outgroup(node, outgroup_branch_length=original_length * fraction)
    tree.root.name = 'Outgroup-supported ingroup ancestor'
    return tree, {
        'method': 'external_outgroup_attachment_projected_to_original_nj',
        'outgroup_id': outgroup_id,
        'ingroup_split_sizes': sorted(map(len, sides)),
        'original_root_edge_length': original_length,
        'fraction_from_matched_side': fraction,
        'dates_used_to_choose_root': False,
        'ingroup_tip_ids_unchanged': True,
        'ingroup_pairwise_branch_distances_unchanged': True,
    }
