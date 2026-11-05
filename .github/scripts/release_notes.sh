#!/bin/bash
# Print a version's CHANGELOG.md section as GitHub release notes.
# Usage: release_notes.sh <version> [changelog]
#
# The section runs from "## <version> ..." (or "## [<version>] ...", Keep a Changelog
# style) to the next "## " heading, without the heading itself.  Hard-wrapped lines are joined, because GitHub shows each line
# break in a release description as a break on the page.  If the section is
# missing or empty, warns on stderr and prints a pointer to CHANGELOG.md instead,
# so a forgotten entry never blocks a release.
set -euo pipefail

if [ $# -lt 1 ] || [ $# -gt 2 ]; then
    echo "usage: $0 <version> [changelog]" >&2
    exit 2
fi
version="$1"
changelog="${2:-CHANGELOG.md}"

notes="$(awk -v version="$version" '
    function flush() { if (buf != "") print buf; buf = "" }
    /^## / { if (in_section) exit; v = $2; gsub(/^\[|\]$/, "", v); in_section = (v == version); next }
    !in_section { next }
    /^[[:space:]]*$/ { flush(); if (printed) blank = 1; next }
    {
        starts_block = ($0 ~ /^#/ || $0 ~ /^[[:space:]]*([-*+]|[0-9]+\.)[[:space:]]/)
        if (buf != "" && !starts_block) {
            line = $0
            sub(/^[[:space:]]+/, "", line)
            buf = buf " " line
            next
        }
        flush()
        if (blank) { print ""; blank = 0 }
        buf = $0
        printed = 1
    }
    END { flush() }
' "$changelog")"

if [ -z "$notes" ]; then
    echo "::warning::no \"## [$version]\" section with any content in $changelog; using boilerplate release notes" >&2
    if [ -n "${GITHUB_REPOSITORY:-}" ] && [ -n "${GITHUB_REF_NAME:-}" ]; then
        notes="See details in [CHANGELOG.md](https://github.com/$GITHUB_REPOSITORY/blob/$GITHUB_REF_NAME/CHANGELOG.md)."
    else
        notes="See details in CHANGELOG.md."
    fi
fi
printf '%s\n' "$notes"
