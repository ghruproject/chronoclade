from copy import deepcopy

from chronoclade import pasteur_profile_normalization as m
from chronoclade.progress import configure_progress
from chronoclade.query_typing import _assignment_result


def test_assignment_summary_uses_actual_parser_statuses_and_does_not_modify_rows(capsys):
    config={'kind':'plincer','scheme':'scgMLST629_S','scheme_version':'fixture'}
    raw=[{'LINcode':['0']*10,'cgST':'1','matches':[]},
         {'LINcode':['0']*6+['*']*4,'cgST':'unassigned','matches':[]},
         {'LINcode':[],'cgST':'unassigned','matches':[]}]
    rows=[_assignment_result(value,config) for value in raw]
    before=deepcopy(rows)
    configure_progress()
    m._assignment_progress(rows)
    assert capsys.readouterr().err=='[typing] cgLIN assignment complete: 1 resolved, 1 provisional, 1 unassigned\n'
    assert rows==before


def test_assignment_summary_respects_quiet(capsys):
    configure_progress(quiet=True)
    rows=[{'sample_id':'PRIVATE_SAMPLE','cglin_status':'complete'}]
    before=deepcopy(rows)
    m._assignment_progress(rows)
    assert capsys.readouterr().err=='' and rows==before
