from pathlib import Path
import sys
import pathlib

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent / "src"))

from unisul.__main__ import main  # ← nicht "src.unisul", sondern Paket "unisul"

if __name__ == "__main__":
    main()
