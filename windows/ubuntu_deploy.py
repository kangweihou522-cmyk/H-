"""Optional Ubuntu deployment using the same versioned D-drive layout as Windows."""
import argparse
import os
from pathlib import Path
import shutil
import shlex
import subprocess
import sys

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from windows.deploy import prepare, activate, read_config, ensure_idle
from windows.station_runtime import stop
from windows.ubuntu_shortcut import create


def install(source, root, port=None):
    source, root = Path(source).resolve(), Path(root).resolve()
    old = read_config(root)
    port = port if port is not None else old.get('port', 8767)
    if not 1024 <= port <= 65535:
        raise ValueError('端口必须为1024–65535。')
    ensure_idle(root)
    if not shutil.which('docker'):
        raise ValueError('Ubuntu 找不到 Docker；请打开 Docker Desktop 并启用此 Ubuntu 的 WSL Integration。')
    subprocess.run(['docker', 'info'], check=True, stdout=subprocess.DEVNULL)
    subprocess.run(['bash', str(source/'windows/wsl_entry.sh'), 'prerequisites', str(root), str(source)], check=True)
    if old:
        stop(root, old['port'])
    windows_root = subprocess.check_output(['wslpath', '-w', str(root)], text=True).strip() if shutil.which('wslpath') else str(root)
    distro = os.environ.get('WSL_DISTRO_NAME') or old.get('distro', 'Ubuntu-24.04')
    proposal = prepare(source, root, windows_root, distro, port)
    cfg = proposal['config']
    subprocess.run(['bash', cfg['program_wsl']+'/windows/wsl_entry.sh', 'install', cfg['root_wsl'],
                    cfg['program_wsl'], cfg['venv_wsl'], str(cfg['port']), cfg['root_windows']], check=True)
    cfg = activate(root, proposal['plan_file'])
    entry = create(root)
    if entry:
        print('以后启动：bash '+shlex.quote(str(entry)), flush=True)
    return cfg


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', default='/mnt/d/转录组分析工作站')
    parser.add_argument('--port', type=int)
    args = parser.parse_args()
    try:
        if args.root.startswith('/mnt/d/') and not Path('/mnt/d').is_mount():
            raise ValueError('Ubuntu 尚未挂载 D 盘（/mnt/d），请先检查 WSL 的磁盘挂载。')
        cfg = install(Path(__file__).resolve().parent.parent, args.root, args.port)
        print('部署完成，正在启动后台。Ubuntu快捷入口见上方输出。', flush=True)
        os.execvp('bash', ['bash', cfg['program_wsl']+'/windows/wsl_entry.sh', 'start', cfg['root_wsl'],
                          cfg['program_wsl'], cfg['venv_wsl'], str(cfg['port']), cfg['root_windows']])
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.exit(1, 'Ubuntu 部署失败：'+str(error)+'\n')
