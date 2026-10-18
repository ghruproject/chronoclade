# Follow-up fixes after independent review

The independent review in COMMAND_REVIEW.md and receipt.json was written before these edits.

Changed chronoclade/tree_stage/cli.py to print the first failed selection reason (including its audit/log path), with Rich markup disabled; failed ensembles still exit 2. Changed chronoclade/execution/ensemble.py to pass completed-run count to the output guide. Changed chronoclade/stage_readme.py to omit the time command and give recovery instructions when a tree ensemble has zero completed runs. Successful tree guides retain their time handoff.

Focused regression checks in tests/test_ensemble_cli_status.py and tests/test_stage_ensembles.py cover actionable literal-path failure output, unchanged failure exit status, zero-completed recovery guidance, and successful-tree handoff preservation.

Validation:

```sh
pixi run pytest tests/test_ensemble_cli_status.py tests/test_stage_ensembles.py -q
```

17 passed in 25.01s.

```sh
pixi run ruff check chronoclade/tree_stage/cli.py chronoclade/execution/ensemble.py chronoclade/stage_readme.py tests/test_ensemble_cli_status.py tests/test_stage_ensembles.py
```

All checks passed.

Manual acceptance rerun:

```sh
pixi run chronoclade tree tmp/beta-first-user-review/run/cgmlst/cglin_7bbf5faa1480/selection_bundle/ensemble.json --out tmp/beta-first-user-review/run/tree-fixed --threads 4 --randomisations 100 --ensemble-jobs 1
```

Exit 2. Stdout includes `Failed selection-001: Selected assembly acquisition failed; no replacements chosen. See .../assembly_audit.json`. The new tree-fixed/README.md says `No tree selections completed. Fix the failures before running time.` and gives recovery guidance without a time command.

No Git mutations or scheduler submissions were performed.
