from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))  # ← src zum Importpfad hinzufügen

from hmsss.__main__ import main  # ← nicht "src.hmsss", sondern Paket "hmsss"

if __name__ == "__main__":
    main()
