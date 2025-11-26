#!/bin/bash
set -euo pipefail

GMX=${GMX:-gmx}
TRAJ=${1:-../md/md.xtc}
TOP=${2:-../md/md.tpr}
OUTDIR=${3:-../analysis}
mkdir -p "$OUTDIR"

$GMX cluster -s "$TOP" -f "$TRAJ" -method gromos -cutoff 0.25 -o "$OUTDIR/cluster.xpm" -g "$OUTDIR/cluster.log"
