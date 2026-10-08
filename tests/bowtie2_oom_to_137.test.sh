#!/usr/bin/env bash
# Checks bin/bowtie2_oom_to_137 (#107) with fake commands that print what
# bowtie2's wrapper and MetaPhlAn print. Run: bash tests/bowtie2_oom_to_137.test.sh
set -u
wrap="$(cd "$(dirname "$0")/.." && pwd)/bin/bowtie2_oom_to_137"
cd "$(mktemp -d)"
fail=0

check() { # name, expected exit, command...
    local name=$1 want=$2; shift 2
    "$wrap" "$@" >out.txt 2>err.txt
    local got=$?
    if [ "$got" -eq "$want" ]; then echo "ok   $name"; else echo "FAIL $name: exit $got, want $want"; fail=1; fi
}

check "OOM via shell (exited with value 137) -> 137" 137 \
    bash -c 'echo "(ERR): bowtie2-align exited with value 137" >&2; echo "[Error] Error while running bowtie2." >&2; exit 1'
check "OOM via signal (died with signal 9) -> 137" 137 \
    bash -c 'echo "(ERR): bowtie2-align died with signal 9 (KILL) " >&2; exit 1'
check "other bowtie2 failure keeps its exit" 1 \
    bash -c 'echo "(ERR): bowtie2-align exited with value 1" >&2; exit 1'
check "signal 15 is not an OOM" 1 \
    bash -c 'echo "(ERR): bowtie2-align died with signal 15 (TERM) " >&2; exit 1'
check "unrelated failure keeps its exit" 2 \
    bash -c 'echo "some other error" >&2; exit 2'
check "success is untouched" 0 \
    bash -c 'echo profile; echo "(ERR): bowtie2-align exited with value 137" >&2; exit 0'
grep -qx profile out.txt || { echo "FAIL stdout not passed through"; fail=1; }
grep -q "exited with value 137" err.txt || { echo "FAIL stderr not passed through"; fail=1; }
[ -z "$(ls -A | grep bowtie2_oom)" ] || { echo "FAIL temp file left behind"; fail=1; }

exit $fail
