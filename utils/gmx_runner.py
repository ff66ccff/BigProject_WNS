import os
import subprocess

def run_gmx_mdrun_safe(deffnm, gmx_cmd="gmx"):
    """
    Runs gmx mdrun with automatic checkpoint detection.
    """
    cpt_file = f"{deffnm}.cpt"
    cmd = [gmx_cmd, "mdrun", "-deffnm", deffnm]

    if os.path.exists(cpt_file):
        print(f"🔄 Found checkpoint '{cpt_file}'. Resuming simulation...")
        # -cpi: Continue Previous simulation
        # -append: Append output to log/edr/trr files instead of creating new ones
        cmd.extend(["-cpi", cpt_file, "-append"])
    else:
        print(f"🚀 Starting new simulation for '{deffnm}'...")

    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"❌ MD Simulation failed for {deffnm}")
        raise e