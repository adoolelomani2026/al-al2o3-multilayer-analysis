from __future__ import annotations

import subprocess
import sys
from pathlib import Path


SCRIPTS = [
    "analyze_stem_eels.py",
    "analyze_uvvis.py",
    "analyze_cl.py",
    "analyze_pl.py",
]


def main() -> None:
    scripts_dir = Path(__file__).resolve().parent
    for script in SCRIPTS:
        print(f"Running {script}")
        subprocess.run([sys.executable, str(scripts_dir / script)], check=True)


if __name__ == "__main__":
    main()
