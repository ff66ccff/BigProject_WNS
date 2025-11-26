param (
    [string]$MGLToolsRoot = "D:\MGLTools",
    [string]$InputPdb = "..\pdb\protein_clean.pdb",
    [string]$OutputPdbqt = "..\pdbqt\protein.pdbqt",
    [int]$Seed = 20231129
)

$prepareScript = Join-Path $MGLToolsRoot "MGLTools-1.5.7\Lib\site-packages\AutoDockTools\Utilities24\prepare_receptor4.py"
$pythonExe = Join-Path $MGLToolsRoot "MGLTools-1.5.7\mglpython.exe"

if (-not (Test-Path $pythonExe)) {
    throw "Unable to find mglpython.exe at $pythonExe. Adjust the MGLToolsRoot parameter."
}

if (-not (Test-Path $prepareScript)) {
    throw "prepare_receptor4.py not found at $prepareScript. Verify the MGLTools installation."
}

if (-not (Test-Path $InputPdb)) {
    throw "Input PDB not found at $InputPdb. Run preprocess_pdb.py first."
}

$cmd = @(
    "-s", $prepareScript,
    "-r", (Resolve-Path $InputPdb),
    "-o", (Resolve-Path $OutputPdbqt),
    "-A", "checkhydrogens",
    "-U", "nphs_lps",
    "-e", "seed=$Seed"
)

& $pythonExe @cmd
if ($LASTEXITCODE -ne 0) {
    throw "prepare_receptor4.py exited with code $LASTEXITCODE"
}
