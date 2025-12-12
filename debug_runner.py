
import sys
import traceback
from pathlib import Path

# Add scripts directory to path to import modules
sys.path.append(str(Path("scripts").resolve()))

# Import the main function
from run_full_wrap_n_shake import main

if __name__ == "__main__":
    try:
        # Run with arguments: --config scripts/config.yml --reset
        sys.argv = ["run_full_wrap_n_shake.py", "--config", "scripts/config.yml", "--reset"]
        main()
    except Exception:
        with open("debug_error.log", "w", encoding="utf-8") as f:
            f.write("\n\n=== CAPTURED EXCEPTION ===\n")
            traceback.print_exc(file=f)
        print("Exception captured to debug_error.log")
