# Releasing lightsim2grid

A release is a pull request from `dev_x.y.z` into `master`. The mechanical parts are in
`utils/release.py`; the parts no script can check are a checklist in the pull request
description. Below, `x.y.z` is the version being released and `x.y.z'` the next one
(usually `z + 1`).

## Before merging

1. On `dev_x.y.z`:

   ```
   python3 utils/release.py prepare x.y.z            # today's date, or --date YYYY-MM-DD
   ```

   It writes the final version into `lightsim2grid/__init__.py` and `docs/conf.py`, dates the
   `[x.y.z]` heading of `CHANGELOG.rst` and sorts its entries by tag
   (BREAKING > DEPRECATED > FIXED > ADDED > IMPROVED > UPDATED, an entry keeping its place
   among those of its tag). Commit and push.

2. Open the pull request into `master` with the release template:
   `https://github.com/Grid2op/lightsim2grid/compare/master...dev_x.y.z?template=release.md`.
   If it was opened without it, the `Release checklist` workflow appends the checklist to
   the description.

3. Go through the checklist. Two workflows run on every pull request into `master`, and
   their three checks fail until:

   - **Version and changelog** (`release_check.yml`): `python3 utils/release.py check` passes;
   - **Docs build without warnings** (`release_check.yml`): Sphinx, installed as Read the
     Docs installs it, with `-W`;
   - **Release checklist** (`release_checklist.yml`): every box ticked, *Release notes* and
     *Announcement* written. Ticking a box reruns it.

   The boxes:

   - benchmarks run (`benchmarks/run_all_benchmarks.sh`) and `docs/benchmarks*.rst` updated;
   - docs updated for the user-facing changes of the changelog;
   - the Read the Docs build of the pull request read;
   - the full description of the changes (it becomes the GitHub release text);
   - the short Discord / Slack announcement.

4. The usual CI (`main.yml`, the C++ jobs) green as well.

## After merging

5. **Tag.** On the merge commit of `master`:

   ```
   git fetch origin master
   git tag -a vx.y.z origin/master -m "lightsim2grid x.y.z"
   git push origin vx.y.z
   ```

   The tag push runs `main.yml`, which builds every wheel and the sdist into the
   `lightsim2grid-wheels` artifact. Wait for it to be green.

6. **GitHub release.** Create a release on the tag `vx.y.z`, titled `x.y.z`, with the
   *Release notes* of the pull request as its text, and publish it. Publishing runs
   `publish.yml`, which checks the tag against `__version__`, re-counts the wheels and
   uploads them to PyPI.

7. **Watch PyPI.** Follow the `Publish to PyPI` run, then check that
   https://pypi.org/project/lightsim2grid/x.y.z/ lists every wheel and that
   `pip install lightsim2grid==x.y.z` works in a fresh environment. If the run failed, read
   why before rerunning it (see the comments in `publish.yml`).

8. **Announce.** Post the *Announcement* of the pull request on Discord and Slack, once PyPI
   has the version.

9. **Open the next development branch.**

   ```
   git checkout -b dev_x.y.z' origin/master
   python3 utils/release.py start-dev x.y.z'
   git commit -as -m "Start x.y.z' development"
   git push -u origin dev_x.y.z'
   ```

   `start-dev` writes `x.y.z'.dev0` into `lightsim2grid/__init__.py` and `docs/conf.py`,
   opens an empty `[x.y.z'] 20xx-xx-xx` section in `CHANGELOG.rst`, and makes `dev_x.y.z'`
   the default base of `utils/pr_diff_stats.py`. Make `dev_x.y.z'` the default branch of
   the repository on GitHub if it is the one pull requests should target.

10. **Point the services at the new branch.** None of these is in the repository; each is
    a setting of the service:

    - **Read the Docs**: *Versions* → activate `dev_x.y.z'` (and check that `stable` now
      points to `vx.y.z`);
    - **SonarQube**: add `dev_x.y.z'` as a long-lived branch of the project;
    - **Codacy**: *Settings → Branches* → enable `dev_x.y.z'`.

## The `master` branch protection

The checks above only block a merge if `master` requires them. In the repository's
*Settings → Rules → Rulesets* (or *Settings → Branches*), for the branch `master`:

- require a pull request before merging;
- require status checks to pass, with `Version and changelog`,
  `Docs build without warnings` and `Release checklist` (plus the `main.yml` jobs wanted
  green). GitHub only offers a check once it has run, so open the release pull request
  first;
- block force pushes and deletion.

Keep the checklist in its own workflow: a job skipped on a description edit would count
as passed and hide a failing check of the same commit.
