"""Build a wheel and exercise it in a new environment, on any CI host."""
import os
from pathlib import Path
import subprocess
import sys
import venv

repo = Path(__file__).resolve().parents[1]
root = Path(os.environ['RUNNER_TEMP'])/'gcf-installed-report'
root.mkdir()
subprocess.run([sys.executable, '-m', 'pip', 'wheel', str(repo), '--no-deps',
                '--wheel-dir', str(root/'wheels')], check=True)
venv.EnvBuilder(with_pip=True).create(root/'env')
python = root/'env'/('Scripts/python.exe' if os.name == 'nt' else 'bin/python')
wheel, = (root/'wheels').glob('*.whl')
subprocess.run([str(python), '-m', 'pip', 'install', '--no-index', '--no-deps', str(wheel)], check=True)
subprocess.run([str(python), str(repo/'tests/installed_report_smoke.py'), '--out',
                str(root/'proof'), '--runtime', str(root.parent/'gcf-runtime')], cwd=root, check=True)
