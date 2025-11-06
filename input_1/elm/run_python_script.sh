#!/bin/bash

WORKDIR=$PWD
SCRIPT=""

USER_ARGS=()  # will hold args for python after `--`
FORWARD=0


while [[ "$#" -gt 0 ]]; do
    if [[ $FORWARD -eq 1 ]]; then
        # after `--`, just collect everything verbatim
        USER_ARGS+=("$1")
        shift
        continue
    fi

    case "$1" in
        --workdir)
            WORKDIR="$2"
            shift 2
            ;;
        --script)
            SCRIPT="$2"
            shift 2
            ;;
        --)
            # everything after this goes to python, not to sbatch wrapper
            FORWARD=1
            shift
            ;;
        *)
            echo "Unknown parameter to run_python.slurm: $1"
            echo "Did you mean to put arguments to python after -- ?"
            exit 1
            ;;
    esac
done

# --- Sanity checks ---
if [[ -z "$SCRIPT" ]]; then
    echo "ERROR: you must specify --script <file.py>"
    exit 1
fi

if [[ ! -d "$WORKDIR" ]]; then
    echo "ERROR: workdir '$WORKDIR' does not exist"
    exit 1
fi

# --- Activate environment (adapt if your cluster uses module load etc.) ---
source ~/.bashrc
conda activate biotools

# --- Move to working directory ---
cd "$WORKDIR" || exit 1
echo "[INFO] Workdir: $(pwd)"
echo "[INFO] Script:  $SCRIPT"
echo "[INFO] Python:  $(which python)"
echo "[INFO] Passing to python: ${USER_ARGS[@]}"

# --- Run the python script with forwarded args ---
python "$SCRIPT" "${USER_ARGS[@]}"
