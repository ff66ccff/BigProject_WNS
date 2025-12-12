#!/usr/bin/env python3
import sys
sys.path.append('scripts')
from wrap_n_shake_docking import extract_best_pose
from pathlib import Path

try:
    extract_best_pose(Path('autodock_runs/wrapper_101.dlg'), Path('autodock_runs/test_extracted.pdbqt'))
    print('Extraction successful')
except Exception as e:
    print(f'Error: {e}')