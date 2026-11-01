<!-- release-checklist -->
<!--
Release pull request: dev_x.y.z -> master. The whole procedure is in RELEASE.md.
The "release checklist" CI job fails until every box below is ticked and the two
text sections are written; it reruns each time this description is edited.
-->

<!-- Diff breakdown: python3 utils/pr_diff_stats.py origin/master -->

## Checked by CI

- `python3 utils/release.py check`: final `x.y.z` version in `lightsim2grid/__init__.py` and
  `docs/conf.py`, matching the top `CHANGELOG.rst` section, dated, and sorted
  BREAKING > DEPRECATED > FIXED > ADDED > IMPROVED > UPDATED.
  `python3 utils/release.py prepare x.y.z` does all of it.
- The documentation builds with no Sphinx warning (`release_check.yml`, job `docs`).

## Before merging

- [ ] Benchmarks run (`benchmarks/run_all_benchmarks.sh`) and the benchmark pages of the docs (`docs/benchmarks*.rst`) updated
- [ ] Docs updated for every user-facing change of the changelog
- [ ] Read the Docs build of this pull request opened and read (rendering, not only warnings)
- [ ] Release notes written below: the full description of the changes
- [ ] Discord / Slack announcement written below

## Release notes

<!--
The full description of the changes. It becomes the text of the GitHub release:
the highlights first, then what breaks and how to migrate, then a link to the changelog.
-->

## Announcement (Discord / Slack)

<!-- A few lines, a link to the release and to the docs. Posted once PyPI has the version. -->
