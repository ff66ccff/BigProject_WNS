#!/bin/bash
set -euo pipefail

GMX=${GMX:-gmx}
# INPUT is relative to the current directory (working_dir), not to WORKDIR
INPUT=${1:-complex/complex_filtered.pdb}
WORKDIR=${2:-gmx}
LIGAND_MOL2=${3:-data/intermediate1.mol2}  # Ligand file for parameterization

# Resolve INPUT to absolute path before changing directory
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CURRENT_DIR="$(pwd)"

# If INPUT is relative, make it absolute relative to current directory
if [[ ! "$INPUT" = /* ]]; then
    INPUT="$CURRENT_DIR/$INPUT"
fi

# If LIGAND_MOL2 is relative, make it absolute
if [[ ! "$LIGAND_MOL2" = /* ]]; then
    LIGAND_MOL2="$CURRENT_DIR/$LIGAND_MOL2"
fi

mkdir -p "$WORKDIR"
pushd "$WORKDIR" >/dev/null

# ============================================================
# Step 0: Extract protein-only PDB and parameterize ligand
# ============================================================
echo ">>> Extracting protein from complex PDB..."
PROTEIN_PDB="protein_only.pdb"
grep -E "^(ATOM|TER)" "$INPUT" | grep -v "UNL\|LIG\|MOL" > "$PROTEIN_PDB" || true

# Check if ligand exists in complex
if grep -qE "UNL|LIG|MOL" "$INPUT"; then
    echo ">>> Ligand detected in complex, generating parameters with ACPYPE..."
    
    # Activate conda and run ACPYPE
    source ~/miniconda3/etc/profile.d/conda.sh
    conda activate amber
    
    # Run ACPYPE on ligand mol2 file with -f to force processing even with close atoms
    LIGAND_NAME="LIG"
    echo ">>> Running ACPYPE on $LIGAND_MOL2..."
    if ! acpype -i "$LIGAND_MOL2" -n 0 -a gaff2 -o gmx -f 2>&1; then
        echo ">>> Retrying ACPYPE with auto charge detection..."
        if ! acpype -i "$LIGAND_MOL2" -a gaff2 -o gmx -f 2>&1; then
            echo ">>> Retrying ACPYPE with gaff force field..."
            acpype -i "$LIGAND_MOL2" -a gaff -o gmx -f 2>&1 || {
                echo "ERROR: All ACPYPE attempts failed!"
                HAS_LIGAND=0
            }
        fi
    fi
    
    # Find ACPYPE output directory
    ACPYPE_DIR=$(ls -d *.acpype 2>/dev/null | head -1)
    if [[ -n "$ACPYPE_DIR" && -d "$ACPYPE_DIR" ]]; then
        echo ">>> ACPYPE output found in $ACPYPE_DIR"
        cp "$ACPYPE_DIR"/*_GMX.itp ./ligand.itp 2>/dev/null || cp "$ACPYPE_DIR"/*.itp ./ligand.itp
        cp "$ACPYPE_DIR"/*_GMX.gro ./ligand.gro 2>/dev/null || cp "$ACPYPE_DIR"/*.gro ./ligand.gro
        # Get ligand residue name from ITP
        LIGAND_RESNAME=$(grep "^\[ moleculetype \]" -A2 ligand.itp | tail -1 | awk '{print $1}')
        echo ">>> Ligand residue name: $LIGAND_RESNAME"
        HAS_LIGAND=1
    else
        echo "Warning: ACPYPE output not found, proceeding with protein only"
        HAS_LIGAND=0
    fi
else
    echo ">>> No ligand detected in complex, proceeding with protein only"
    HAS_LIGAND=0
fi

# ============================================================
# Step 1: Generate protein topology
# ============================================================
echo ">>> Running pdb2gmx on protein..."
$GMX pdb2gmx -f "$PROTEIN_PDB" -o protein.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh

# ============================================================
# Step 2: Combine protein and ligand if ligand exists
# ============================================================
if [[ "${HAS_LIGAND:-0}" -eq 1 ]]; then
    echo ">>> Combining protein and ligand..."
    
    # Get protein atom count (exclude last 2 lines: box vectors line and blank)
    PROTEIN_ATOMS=$(head -2 protein.gro | tail -1)
    LIGAND_ATOMS=$(head -2 ligand.gro | tail -1)
    TOTAL_ATOMS=$((PROTEIN_ATOMS + LIGAND_ATOMS))
    
    # Create combined GRO file
    head -1 protein.gro > complex.gro
    echo "$TOTAL_ATOMS" >> complex.gro
    # Add protein coordinates (skip first 2 lines: title and atom count, skip last line: box)
    tail -n +3 protein.gro | head -n -1 >> complex.gro
    # Add ligand coordinates (skip first 2 lines, skip last line)
    tail -n +3 ligand.gro | head -n -1 >> complex.gro
    # Add box vectors from protein
    tail -1 protein.gro >> complex.gro
    
    # Modify topology to include ligand
    # Insert ligand ITP include after the force field include
    sed -i '/^#include.*forcefield\.itp/a #include "ligand.itp"' topol.top
    # Add ligand to [ molecules ] section
    echo "$LIGAND_RESNAME    1" >> topol.top
    
    echo ">>> Combined complex.gro created with $TOTAL_ATOMS atoms"
else
    # Just use protein as complex
    cp protein.gro complex.gro
fi

# ============================================================
# Step 3: System preparation (box, solvation, ions)
# ============================================================
echo ">>> Setting up simulation box..."
$GMX editconf -f complex.gro -o boxed.gro -c -d 1.0 -bt cubic

echo ">>> Solvating system..."
$GMX solvate -cp boxed.gro -cs tip3p.gro -o solvated.gro -p topol.top

echo ">>> Adding ions..."
$GMX grompp -f ../mdp/ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1
$GMX genion -s ions.tpr -o ionized.gro -neutral -conc 0.15 <<<'SOL'

# ============================================================
# Step 4: Energy minimization
# ============================================================
echo ">>> Running energy minimization..."
$GMX grompp -f ../mdp/em.mdp -c ionized.gro -p topol.top -o em.tpr -maxwarn 1
$GMX mdrun -deffnm em

# ============================================================
# Step 5: Equilibration (NVT and NPT)
# ============================================================
echo ">>> Running NVT equilibration..."
$GMX grompp -f ../mdp/nvt.mdp -c em.gro -p topol.top -o nvt.tpr -maxwarn 1
$GMX mdrun -deffnm nvt

echo ">>> Running NPT equilibration..."
$GMX grompp -f ../mdp/npt.mdp -c nvt.gro -p topol.top -o npt.tpr -maxwarn 1
$GMX mdrun -deffnm npt

# ============================================================
# Step 6: Production MD
# ============================================================
echo ">>> Running production MD..."
$GMX grompp -f ../mdp/md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 1
$GMX mdrun -deffnm md

echo ">>> GROMACS pipeline completed successfully!"
popd >/dev/null
