#!/bin/bash
set -euo pipefail

INPUT=${1:-../data/ligand.mol2}
BASENAME=${2:-ligand}

# Detect input format from extension
case "${INPUT##*.}" in
    pdb|PDB)
        INPUT_FORMAT="pdb"
        ;;
    mol2|MOL2)
        INPUT_FORMAT="mol2"
        ;;
    *)
        INPUT_FORMAT="mol2"
        ;;
esac

if [ -n "${AMBERHOME:-}" ] && [ -f "$AMBERHOME/amber.sh" ]; then
    # shellcheck disable=SC1090
    source "$AMBERHOME/amber.sh"
elif [ -f "/usr/local/amber/amber.sh" ]; then
    # shellcheck disable=SC1090
    source "/usr/local/amber/amber.sh"
elif [ -f "$HOME/amber/amber.sh" ]; then
    # shellcheck disable=SC1090
    source "$HOME/amber/amber.sh"
fi

command -v antechamber >/dev/null 2>&1 || {
    echo "Error: antechamber not found in PATH. Ensure AmberTools environment is sourced." >&2
    exit 127
}

antechamber -i "$INPUT" -fi "$INPUT_FORMAT" -o "${BASENAME}_gaff.mol2" -fo mol2 -c bcc
parmchk2 -i "${BASENAME}_gaff.mol2" -f mol2 -o "${BASENAME}.frcmod"
acpype -i "${BASENAME}_gaff.mol2" -c bcc
