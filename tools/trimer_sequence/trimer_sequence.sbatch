#!/bin/bash
#SBATCH --job-name=trimer_sequence
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

set -uo pipefail
source "${SAPIA_PRELUDE:?}"

set +u
source "${SAPIA_ACTIVATE_TRIMER_SEQUENCE:?set SAPIA_ACTIVATE_TRIMER_SEQUENCE in your .env}"
set -u

name=$(echo "$SAPIA_LINE" | cut -f1)
seq_path=$(echo "$SAPIA_LINE" | cut -f2)
test -f "$seq_path" && echo "OK: $name" || { echo "ERROR: $seq_path not found"; exit 1; }
