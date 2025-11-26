#!/bin/bash
set -euo pipefail

GMX=${GMX:-gmx}
INPUT=${1:-../complex/complex_filtered.pdb}
WORKDIR=${2:-../gmx}
mkdir -p "$WORKDIR"

pushd "$WORKDIR" >/dev/null

$GMX pdb2gmx -f "$INPUT" -o complex.gro -p topol.top -ff amber99sb-ildn
$GMX editconf -f complex.gro -o boxed.gro -c -d 1.0 -bt cubic
$GMX solvate -cp boxed.gro -cs tip3p.gro -o solvated.gro -p topol.top
$GMX grompp -f ../mdp/ions.mdp -c solvated.gro -p topol.top -o ions.tpr
$GMX genion -s ions.tpr -o ionized.gro -neutral -conc 0.15 <<<'SOL'
$GMX grompp -f ../mdp/em.mdp -c ionized.gro -p topol.top -o em.tpr
$GMX mdrun -deffnm em
$GMX grompp -f ../mdp/nvt.mdp -c em.gro -p topol.top -o nvt.tpr
$GMX mdrun -deffnm nvt
$GMX grompp -f ../mdp/npt.mdp -c nvt.gro -p topol.top -o npt.tpr
$GMX mdrun -deffnm npt
$GMX grompp -f ../mdp/md.mdp -c npt.gro -p topol.top -o md.tpr
$GMX mdrun -deffnm md

popd >/dev/null
