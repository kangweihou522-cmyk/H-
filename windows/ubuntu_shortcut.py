"""Install a direct Ubuntu entry without Windows process interop."""
from pathlib import Path
import json
import shlex
import sys


def create(root, home=None):
    root = Path(root)
    cfg = json.loads((root/'配置/workstation.json').read_text())
    file = Path(home or Path.home())/'启动转录组工作站.sh'
    content = ('#!/usr/bin/env bash\nset -euo pipefail\n'
               f'RNASEQ_STATION={shlex.quote(str(root))}\n'
               'mapfile -t rnaseq_paths < <(python3 -c '\
               "'import sys,json;from pathlib import Path;c=json.loads((Path(sys.argv[1])/\"配置/workstation.json\").read_text(encoding=\"utf-8-sig\"));[print(c[k]) for k in (\"program_wsl\",\"venv_wsl\",\"port\",\"root_windows\")]' "\
               '"$RNASEQ_STATION")\n'
               'exec bash "${rnaseq_paths[0]}/windows/wsl_entry.sh" start "$RNASEQ_STATION" "${rnaseq_paths[0]}" "${rnaseq_paths[1]}" "${rnaseq_paths[2]}" "${rnaseq_paths[3]}"\n')
    old_bridge = ('if command -v powershell.exe >/dev/null && command -v wslpath >/dev/null; then\n'
                  '  exec powershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass -File "$(wslpath -w "$RNASEQ_STATION/station.ps1")" -Action Start\n'
                  'fi\n')
    legacy = content.replace('mapfile -t rnaseq_paths', old_bridge+'mapfile -t rnaseq_paths', 1)
    if file.exists() and file.read_text() not in {content, legacy}:
        # Preserve a user-edited entry or another station, use a distinct new entry.
        file = file.with_name('启动D盘转录组工作站.sh')
        if file.exists() and file.read_text() != content:
            print('Existing Ubuntu entries preserved; use the D-drive Windows launcher.')
            return None
    file.write_text(content, encoding='utf-8');file.chmod(0o755)
    print(str(file))
    return file


if __name__ == '__main__':
    create(sys.argv[1])
