#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# QA LAYER 1 - static text scan
#
# Fails if any plain-text file git would publish contains a participant-level
# clinical column name, a personal identifier pattern, or a credential.
# Column *names* are the signal: a row-per-participant table announces itself
# in its header.
#
# Usage:  bash scripts/qa/scan_pii.sh
# Exit 0 = clean, 1 = finding.
# ---------------------------------------------------------------------------
set -uo pipefail
cd "$(git rev-parse --show-toplevel)" || exit 1

# Header tokens that indicate a per-participant clinical table.
# Only variables the published article itself names (SI Table 6: "~ Group +
# vaginal pH + antibiotic + BV history + age"). Cohort variables that the
# article does not mention are guarded in the private authoring repository's
# copy of this scanner, so that this public file does not enumerate them.
# Word-anchored so 'Age' does not match 'Average'/'page'. Aligned with the
# BANNED list in scan_rda.R; Group_BV is excluded there and here because the
# BV/Control label is published.
BANNED='Vaginal[_ ]?pH|BV[_ ]?History|pH_z|Age_z|Base_Sample_ID|Sample_ID|(^|[^A-Za-z0-9_])(Age|Antibiotic)([^A-Za-z0-9_]|$)|personnummer|date_of_birth'
# Secrets / credentials.
# Accept ':' as well as '=': YAML config files write "API_KEY : value".
SECRETS='BEGIN [A-Z ]*PRIVATE KEY|api[_-]?key[[:space:]]*[:=]|password[[:space:]]*[:=]|secret[[:space:]]*[:=]|token[[:space:]]*[:=][[:space:]]*[A-Za-z0-9]{16,}|ghp_[A-Za-z0-9]{20,}|AKIA[0-9A-Z]{16}'
# Addresses other than the corresponding author's in README.md.
EMAILS='[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}'
# Local absolute paths that would leak machine structure.
# Local machine paths and HPC site paths (the latter leak allocation layout).
LOCALPATH='/Users/[a-zA-Z]|/home/[a-zA-Z]|C:\\\\Users|/crex/|/gorilla/|/proj/[a-z]|/sw/data/'

# Scripts legitimately reference these column names in code; data files must not.
# Markdown documentation names these variables in prose by design (legends,
# README, controlled-access notes). The risk is tabular data, so .md is exempt
# from the column-name check but still subject to the secret/path/email checks.
CODE_EXEMPT='^scripts/.*\.(R|sh)$|\.md$|^\.gitignore$|^\.githooks/|^\.github/'
# The QA tooling and the hook contain these patterns as literals by design, so
# they are exempt from the secret / path / email checks (not from BANNED).
SELF_EXEMPT='^scripts/qa/|^\.githooks/'

fail=0
report() { printf '  %-7s %s\n' "$1" "$2"; fail=1; }

echo "QA L1 - static text scan"
files=$(git ls-files --cached --others --exclude-standard)
n=0

while IFS= read -r f; do
  [ -f "$f" ] || continue
  # skip binaries
  # Binaries that layer 2 does not cover (it handles .rda/.rds/.xlsx) are
  # scanned through `strings` rather than skipped outright.
  if file --mime "$f" 2>/dev/null | grep -q 'charset=binary'; then
    case "$f" in
      *.rda|*.RData|*.rds|*.xlsx) continue ;;
      *) if strings "$f" 2>/dev/null | grep -qiE "$BANNED"; then
           report "BINARY" "$f (clinical name found in binary content)"
         fi
         continue ;;
    esac
  fi
  # Raised from 5 MB: the previous cap silently skipped any larger table.
  [ "$(wc -c <"$f")" -gt 268435456 ] && { report "SKIP" "$f (>256MB, not scanned)"; continue; }
  n=$((n+1))

  if ! printf '%s' "$f" | grep -qE "$CODE_EXEMPT"; then
    hits=$(grep -ioE "$BANNED" "$f" 2>/dev/null | sort -u | tr '\n' ',' | sed 's/,$//')
    [ -n "$hits" ] && report "CLINIC" "$f -> $hits"
  fi
  if ! printf '%s' "$f" | grep -qE "$SELF_EXEMPT"; then
    s=$(grep -ioE "$SECRETS" "$f" 2>/dev/null | head -1)
    [ -n "$s" ] && report "SECRET" "$f"
    p=$(grep -oE "$LOCALPATH" "$f" 2>/dev/null | head -1)
    [ -n "$p" ] && report "PATH" "$f -> $p"
    # README.md carries the corresponding author's address by design.
    if [ "$f" != "README.md" ]; then
      e=$(grep -oE "$EMAILS" "$f" 2>/dev/null | head -1)
      [ -n "$e" ] && report "EMAIL" "$f"
    fi
  fi
done <<< "$files"

echo "scanned $n text file(s) git would publish"
if [ "$fail" -ne 0 ]; then
  echo "FAIL - see findings above."
  exit 1
fi
echo "PASS - no participant identifiers, secrets, or local paths."
exit 0
