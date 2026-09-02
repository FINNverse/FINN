#!/usr/bin/env python3
"""PreToolUse (Bash) hook: block bulk `git add` in the FINN repo.

`git add -A` / `--all` / `.` has previously swept unpublished manuscripts and
personal scratch into commits. This hook refuses those forms and tells Claude to
stage specific paths instead. It fails OPEN (exit 0) on any parsing problem, so a
missing python3 or an unexpected payload never blocks legitimate work — the real
safety net is per-file staging plus CI on PRs.

Exit codes: 0 = allow, 2 = block (stderr is shown to Claude).
"""
import sys, json, re

try:
    data = json.load(sys.stdin)
except Exception:
    sys.exit(0)  # fail open

if data.get("tool_name") not in (None, "Bash"):
    sys.exit(0)

cmd = (data.get("tool_input") or {}).get("command", "") or ""

# Match `git add -A/--all/.` only at a real command position: the start of the
# command, or right after a shell separator (; & | newline ( ). This ignores the
# phrase when it appears as quoted data (echo/grep/commit messages), so only an
# actual bulk-add invocation is blocked.
pattern = re.compile(
    r"(?:^|[;&|\n(])\s*(?:cd\b[^\n;&|]*?(?:&&|;)\s*)?"   # optional leading `cd ... &&`
    r"git\s+add\b[^\n&|;]*?(\s-A\b|\s--all\b|\s\.(?:\s|$))"
)
if pattern.search(cmd):
    sys.stderr.write(
        "Blocked: bulk `git add -A/--all/.` is forbidden in the FINN repo.\n"
        "It has swept unpublished manuscripts and personal scratch into commits.\n"
        "Stage specific paths instead, e.g. `git add R/finn-class.R man/finn.Rd`.\n"
    )
    sys.exit(2)

sys.exit(0)
