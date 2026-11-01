#!/usr/bin/env python3
# Copyright (c) 2026, RTE (https://www.rte-france.com)
# See AUTHORS.txt
# This Source Code Form is subject to the terms of the Mozilla Public License, version 2.0.
# If a copy of the Mozilla Public License, version 2.0 was not distributed with this file,
# you can obtain one at http://mozilla.org/MPL/2.0/.
# SPDX-License-Identifier: MPL-2.0
# This file is part of LightSim2grid, LightSim2grid implements a c++ backend targeting the Grid2Op platform.
"""The mechanical parts of a release. See RELEASE.md for the whole procedure.

    python3 utils/release.py prepare 1.3.0 [--date 2026-11-02]
        Before merging into master: final version in lightsim2grid/__init__.py and
        docs/conf.py, date on the changelog heading, changelog section sorted.

    python3 utils/release.py check
        What the "release check" CI job runs on a pull request into master: fails
        when any of the above is not done.

    python3 utils/release.py sort-changelog
        Only sort the entries of the topmost release section of CHANGELOG.rst.

    python3 utils/release.py check-pr-body < body.md
        What the "release checklist" CI job runs on the body of a pull request into
        master: every box of the checklist ticked, the release notes and the
        announcement written.

    python3 utils/release.py start-dev 1.3.1
        After the release, on the new dev_1.3.1 branch: "1.3.1.dev0" in
        lightsim2grid/__init__.py and docs/conf.py, an empty [1.3.1] changelog section,
        and dev_1.3.1 as the default base of utils/pr_diff_stats.py.

Run from anywhere: the paths are taken relative to the repository root.
"""

import argparse
import datetime
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INIT_PY = os.path.join(ROOT, "lightsim2grid", "__init__.py")
CONF_PY = os.path.join(ROOT, "docs", "conf.py")
CHANGELOG = os.path.join(ROOT, "CHANGELOG.rst")
PR_DIFF_STATS = os.path.join(ROOT, "utils", "pr_diff_stats.py")
PR_TEMPLATE = os.path.join(ROOT, ".github", "PULL_REQUEST_TEMPLATE", "release.md")
# The release template carries this line; a body without it never got the checklist.
CHECKLIST_MARKER = "<!-- release-checklist -->"
# Sections of the template that must not be left empty.
MUST_FILL = ["Release notes", "Announcement (Discord / Slack)"]

# The order of the entries inside a release section, by their first tag.
# An entry carrying several tags ("[BREAKING] [DEPRECATED] ...") is sorted by the first.
TAG_ORDER = ["BREAKING", "DEPRECATED", "FIXED", "ADDED", "IMPROVED", "UPDATED"]

INIT_RE = re.compile(r'^__version__\s*=\s*"(?P<v>[^"]+)"$', re.MULTILINE)
CONF_RELEASE_RE = re.compile(r'^release\s*=\s*["\'](?P<v>[^"\']+)["\']$', re.MULTILINE)
CONF_VERSION_RE = re.compile(r'^version\s*=\s*["\'](?P<v>[^"\']+)["\']$', re.MULTILINE)
# "[1.3.0] 2026-11-02" (released) or "[1.3.0] 2026-xx-xx" (in development)
HEADING_RE = re.compile(r"^\[(?P<v>\d+\.\d+\.\d+)\]\s+(?P<date>\S+)\s*$")
FINAL_VERSION_RE = re.compile(r"^\d+\.\d+\.\d+$")
DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")
ENTRY_TAG_RE = re.compile(r"^- \[(?P<tag>[A-Z]+)\]")


def _read(path):
    with open(path, encoding="utf-8") as f:
        return f.read()


def _write(path, text):
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def _sub_once(regex, text, value, path):
    def repl(m):
        start, end = m.span("v")
        return m.group(0)[:start - m.start()] + value + m.group(0)[end - m.start():]
    new, n = regex.subn(repl, text, count=1)
    if n != 1:
        raise SystemExit(f"{path}: pattern {regex.pattern!r} not found")
    return new


def _major_minor(version):
    return ".".join(version.split(".")[:2])


# --- changelog -------------------------------------------------------------------------

def _top_section(lines):
    """(heading index, first body line, end) of the topmost "[x.y.z] date" section."""
    for i, line in enumerate(lines):
        if HEADING_RE.match(line):
            body = i + 1
            if body < len(lines) and set(lines[body].strip()) == {"-"}:
                body += 1  # the rst underline
            end = body
            while end < len(lines) and not HEADING_RE.match(lines[end]):
                end += 1
            return i, body, end
    raise SystemExit(f"{CHANGELOG}: no '[x.y.z] date' heading found")


def _split_entries(body_lines):
    """(entries, trailing) where each entry is the list of its lines.

    An entry starts with "- " and continues on the indented lines below it. The blank
    lines after the last entry are returned apart so that they stay at the end."""
    trailing = []
    while body_lines and not body_lines[-1].strip():
        trailing.insert(0, body_lines.pop())
    entries = []
    for n, line in enumerate(body_lines):
        if line.startswith("- "):
            entries.append([line])
        elif entries and (line.startswith(" ") or not line.strip()):
            entries[-1].append(line)
        else:
            raise SystemExit(f"{CHANGELOG}: cannot read the line {line!r} of the top "
                             "release section as part of an entry (an entry starts with "
                             "'- ', its continuation lines are indented)")
    return entries, trailing


def _entry_rank(entry):
    m = ENTRY_TAG_RE.match(entry[0])
    if m is None or m.group("tag") not in TAG_ORDER:
        return None
    return TAG_ORDER.index(m.group("tag"))


def changelog_errors(lines):
    errors = []
    _, body, end = _top_section(lines)
    entries, _ = _split_entries(lines[body:end])
    last = -1
    for entry in entries:
        rank = _entry_rank(entry)
        if rank is None:
            errors.append(f"CHANGELOG.rst: entry without a known tag ({', '.join(TAG_ORDER)}): "
                          f"{entry[0].strip()[:80]}")
            continue
        if rank < last:
            errors.append("CHANGELOG.rst: the top release section is not sorted "
                          f"({' > '.join(TAG_ORDER)}); first entry out of order: "
                          f"{entry[0].strip()[:80]}. "
                          "Run `python3 utils/release.py sort-changelog`.")
            break
        last = rank
    return errors


def sort_changelog():
    lines = _read(CHANGELOG).split("\n")
    _, body, end = _top_section(lines)
    entries, trailing = _split_entries(lines[body:end])
    for entry in entries:
        if _entry_rank(entry) is None:
            raise SystemExit(f"{CHANGELOG}: entry without a known tag, sort it by hand: "
                             f"{entry[0].strip()[:80]}")
    # sorted() is stable: the entries of one tag keep their order.
    entries = sorted(entries, key=_entry_rank)
    lines[body:end] = [line for entry in entries for line in entry] + trailing
    _write(CHANGELOG, "\n".join(lines))


# --- versions --------------------------------------------------------------------------

def read_versions():
    init = INIT_RE.search(_read(INIT_PY))
    conf = _read(CONF_PY)
    release = CONF_RELEASE_RE.search(conf)
    version = CONF_VERSION_RE.search(conf)
    lines = _read(CHANGELOG).split("\n")
    heading = HEADING_RE.match(lines[_top_section(lines)[0]])
    return {
        "init": init.group("v") if init else None,
        "conf_release": release.group("v") if release else None,
        "conf_version": version.group("v") if version else None,
        "changelog": heading.group("v"),
        "changelog_date": heading.group("date"),
    }


def check():
    v = read_versions()
    errors = []
    pkg = v["init"]
    if pkg is None:
        errors.append("lightsim2grid/__init__.py: no __version__ found")
    elif not FINAL_VERSION_RE.match(pkg):
        errors.append(f"lightsim2grid/__init__.py: __version__ = {pkg!r} is not a release "
                      "version (x.y.z, no .devN / rc suffix)")
    if v["conf_release"] != pkg:
        errors.append(f"docs/conf.py: release = {v['conf_release']!r}, expected {pkg!r}")
    if pkg and v["conf_version"] != _major_minor(pkg):
        errors.append(f"docs/conf.py: version = {v['conf_version']!r}, "
                      f"expected {_major_minor(pkg)!r}")
    if v["changelog"] != pkg_base(pkg):
        errors.append(f"CHANGELOG.rst: the top release section is [{v['changelog']}], "
                      f"expected [{pkg_base(pkg)}]")
    if not DATE_RE.match(v["changelog_date"]):
        errors.append(f"CHANGELOG.rst: the [{v['changelog']}] heading has no date "
                      f"(reads {v['changelog_date']!r}, expected YYYY-MM-DD)")
    errors += changelog_errors(_read(CHANGELOG).split("\n"))

    if errors:
        for e in errors:
            # "::error::" makes GitHub Actions show it as an annotation; harmless locally.
            print(f"::error::{e}" if os.environ.get("GITHUB_ACTIONS") else f"error: {e}")
        print(f"\n{len(errors)} problem(s). `python3 utils/release.py prepare {pkg_base(pkg)}` "
              "fixes the version, the date and the order.")
        return 1
    print(f"ok: version {pkg}, released {v['changelog_date']}, changelog sorted")
    return 0


def pkg_base(version):
    m = re.match(r"^\d+\.\d+\.\d+", version or "")
    return m.group(0) if m else "x.y.z"


def _set_versions(full, short):
    init = _read(INIT_PY)
    _write(INIT_PY, _sub_once(INIT_RE, init, full, INIT_PY))
    conf = _read(CONF_PY)
    conf = _sub_once(CONF_RELEASE_RE, conf, full, CONF_PY)
    conf = _sub_once(CONF_VERSION_RE, conf, short, CONF_PY)
    _write(CONF_PY, conf)


def prepare(version, date):
    if not FINAL_VERSION_RE.match(version):
        raise SystemExit(f"{version!r} is not x.y.z")
    lines = _read(CHANGELOG).split("\n")
    head = _top_section(lines)[0]
    m = HEADING_RE.match(lines[head])
    if m.group("v") != version:
        raise SystemExit(f"CHANGELOG.rst: the top release section is [{m.group('v')}], "
                         f"not [{version}]: fix it by hand first")
    lines[head] = f"[{version}] {date}"
    _write(CHANGELOG, "\n".join(lines))
    sort_changelog()
    _set_versions(version, _major_minor(version))
    return check()


# --- pull request body -----------------------------------------------------------------

def _section_text(body, title):
    """Text under "## title", up to the next "## " heading, html comments removed."""
    m = re.search(r"^##\s+" + re.escape(title) + r"\s*$(?P<t>.*?)(?=^##\s|\Z)",
                  body, re.MULTILINE | re.DOTALL)
    if m is None:
        return None
    return re.sub(r"<!--.*?-->", "", m.group("t"), flags=re.DOTALL).strip()


def pr_body_errors(body):
    body = body.replace("\r\n", "\n")
    if CHECKLIST_MARKER not in body:
        return ["the pull request has no release checklist: its body does not come from "
                ".github/PULL_REQUEST_TEMPLATE/release.md"]
    errors = []
    for line in body.split("\n"):
        if re.match(r"^\s*[-*] \[ \]", line):
            errors.append(f"not ticked: {line.strip()[6:]}")
    for title in MUST_FILL:
        text = _section_text(body, title)
        if text is None:
            errors.append(f"section '## {title}' missing")
        elif not text:
            errors.append(f"section '## {title}' is empty")
    return errors


def check_pr_body(body):
    errors = pr_body_errors(body)
    for e in errors:
        print(f"::error::{e}" if os.environ.get("GITHUB_ACTIONS") else f"error: {e}")
    if not errors:
        print("ok: release checklist complete")
    return 1 if errors else 0


def start_dev(version):
    if not FINAL_VERSION_RE.match(version):
        raise SystemExit(f"{version!r} is not x.y.z")
    _set_versions(f"{version}.dev0", _major_minor(version))

    lines = _read(CHANGELOG).split("\n")
    head = _top_section(lines)[0]
    current = HEADING_RE.match(lines[head]).group("v")
    if current == version:
        print(f"CHANGELOG.rst: [{version}] already there, left as is")
    else:
        title = f"[{version}] 20xx-xx-xx"
        lines[head:head] = [title, "-" * 20, "", ""]
        _write(CHANGELOG, "\n".join(lines))

    stats = _read(PR_DIFF_STATS)
    stats, n = re.subn(r'^DEFAULT_BASE = "[^"]*"$', f'DEFAULT_BASE = "origin/dev_{version}"',
                       stats, count=1, flags=re.MULTILINE)
    if n != 1:
        raise SystemExit(f"{PR_DIFF_STATS}: DEFAULT_BASE not found")
    _write(PR_DIFF_STATS, stats)
    print(f"version {version}.dev0, changelog section [{version}], "
          f"pr_diff_stats.py base origin/dev_{version}")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare", help="final version, dated and sorted changelog")
    p.add_argument("version")
    p.add_argument("--date", default=datetime.date.today().isoformat())
    sub.add_parser("check", help="what CI checks on a pull request into master")
    sub.add_parser("sort-changelog", help="sort the top release section of CHANGELOG.rst")
    sub.add_parser("check-pr-body", help="the release checklist of a PR body (stdin)")
    d = sub.add_parser("start-dev", help="open the next development version")
    d.add_argument("version")
    args = parser.parse_args()

    if args.cmd == "prepare":
        return prepare(args.version, args.date)
    if args.cmd == "check":
        return check()
    if args.cmd == "check-pr-body":
        return check_pr_body(sys.stdin.read())
    if args.cmd == "sort-changelog":
        sort_changelog()
        return 0
    return start_dev(args.version)


if __name__ == "__main__":
    sys.exit(main())
