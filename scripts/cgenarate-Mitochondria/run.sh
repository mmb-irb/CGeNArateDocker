#!/usr/bin/env bash
# Usage: ./run.sh [<FRAMES>] [<STEPS>]
#
# Runs all simulations configured in ./simulations/*.toml files.
#
# If NO_TRACE is set or ovniemu cannot be found, the ovni traces are not
# enabled.
#
# WARN: the frames and steps are only overridden when doing traces.
#
# The behaviour of the program can be tweaked further through environment
# variables:
#
# - TARGET: (default: CGeNArate.exe) name of the make target and binary to run.
# - NO_TRACE: if set, ovni is not used.
# - OVNIEMU: (default: ovniemu-raul) ovniemu program to use.
# - CPUS: cpus to use

set -e

CPUS="${CPUS:-0-55}"
TARGET="${TARGET:-CGeNArate.exe}"
OVNIEMU="${OVNIEMU:-ovniemu-raul}"

make "$TARGET"

doTrace=1
if [ ! -z ${NO_TRACE+x} ]; then
    echo "WARN: disabling trace (NO_TRACE is set)"
    doTrace=0
elif ! command -v "$OVNIEMU" &>/dev/null; then
    echo "WARN: $OVNIEMU not available, disabling trace"
    doTrace=0
fi

if ((doTrace)); then
    export OVNI_TMPDIR="$(mktemp -u /dev/shm/ovni.XXXXXX)"
    export OMP_OVNI=2
    export NOSV_CONFIG_OVERRIDE="instrumentation.version=ovni"

    # We override the frames and steps through this environment variables,
    # so the traces don't get too big.
    export CGEN_FRAMES=${1:-2}
    export CGEN_STEPS=${2:-3}
fi

SIMULATIONS=()

for conf in ./simulations/*.toml; do
    name="$(basename "${conf%.*}")"
    echo "Running $name"

    ((doTrace)) && rm -rf ovni "ovni-$name"
    taskset -c "$CPUS" "./$TARGET" "$conf"
    ((doTrace)) && mv ovni "ovni-$name"

    echo "$name DONE"

    SIMULATIONS+=("$name")
done

python check.py

((doTrace)) || exit 0

echo ALL DONE, running ovniemu

for name in "${SIMULATIONS[@]}"; do
    echo "Simulating $name"
    "$OVNIEMU" "ovni-$name"
done
