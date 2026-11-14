#!/usr/bin/env python3
"""Minimal, dependency-free notebook output stripper for use as a git clean filter.

Configured via:
  git config filter.nbstripout.clean "python3 scripts/git-nbstripout.py"
  git config filter.nbstripout.smudge cat
  git config filter.nbstripout.required true
and in .gitattributes:
  *.ipynb filter=nbstripout

On `git add`/`git commit`, git pipes the notebook's content through this
script's stdin and stores stdout instead. The working-tree copy on disk is
never touched, so outputs are still visible when you open the notebook -
they just never enter git history.

If the input isn't valid notebook JSON (or anything goes wrong), the
original bytes are passed through unchanged so nothing is ever corrupted.
"""
import json
import sys


def strip(nb):
    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        cell["outputs"] = []
        cell["execution_count"] = None
        meta = cell.get("metadata", {})
        for key in ("ExecuteTime", "execution", "collapsed", "scrolled"):
            meta.pop(key, None)
        cell["metadata"] = meta

    # Drop widget state and other execution-only metadata noise that
    # otherwise churns every commit without carrying real information.
    nb.get("metadata", {}).pop("widgets", None)

    return nb


def main():
    raw = sys.stdin.buffer.read()
    try:
        nb = json.loads(raw.decode("utf-8"))
        stripped = strip(nb)
        out = json.dumps(stripped, indent=1, ensure_ascii=False, sort_keys=False)
        sys.stdout.buffer.write(out.encode("utf-8"))
        sys.stdout.buffer.write(b"\n")
    except Exception:
        # Not parseable JSON (or unexpected shape) - pass through untouched.
        sys.stdout.buffer.write(raw)


if __name__ == "__main__":
    main()
