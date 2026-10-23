"""Stage immutable source versions and activate a D-drive workstation after install."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import tempfile
import uuid
from datetime import datetime, timezone

SERVICE = 'mouse-bulk-rnaseq'
DIRECTORIES = ['程序', '原始数据', '参考文件', '分析结果', '导出报告', '日志', '配置']
SOURCE_DIRS = {'scripts', 'static', 'vendor', 'windows', 'docs', 'tests'}
SOURCE_FILES = {'app.py', 'plots.py', 'results.py', 'validation.py', 'requirements.txt',
                'install.sh', 'launch.sh', 'Dockerfile', 'Dockerfile.annotation',
                'README.md', 'VALIDATION.md', 'pytest.ini', '.gitignore', '.gitattributes',
                '启动工作站.bat', '部署到D盘.bat', '检查环境.bat', '更新工作站.bat', '部署到D盘.sh'}
IGNORE = {'.git', '.venv', 'runs', '__pycache__', '.cache', '.pytest_cache'}
REQUIRED = ['app.py', 'install.sh', 'launch.sh', 'requirements.txt',
            'windows/Workstation.ps1', 'windows/Station.ps1', 'windows/wsl_entry.sh']


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_config(root):
    p = root/'配置/workstation.json'
    if not p.exists():
        return {}
    value = json.loads(p.read_text(encoding='utf-8-sig'))
    if value.get('service_id') != SERVICE or value.get('schema') != 1:
        raise ValueError('目标配置不属于本工作站，未覆盖：'+str(p))
    return value


def atomic_write(path, data):
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as f:
        temporary = Path(f.name)
        f.write(data)
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def inventory(source):
    for name in REQUIRED:
        if not (source/name).is_file():
            raise ValueError('源码不完整，请使用新版完整 ZIP：'+name)
    files = {}
    for p in sorted(source.rglob('*')):
        name = p.relative_to(source)
        if any(part in IGNORE for part in name.parts):
            continue
        if name.parts[0] not in SOURCE_DIRS and str(name) not in SOURCE_FILES:
            continue
        if p.is_symlink():
            raise ValueError('源码包含符号链接，未复制：'+str(name))
        if p.is_file():
            files[name.as_posix()] = digest(p.read_bytes())
    identifier = digest(json.dumps(files, sort_keys=True).encode())
    return identifier, files


def ensure_idle(root):
    for p in (root/'分析结果').glob('*/status.json'):
        value = json.loads(p.read_text(encoding='utf-8'))
        if value.get('state') in ('queued', 'running'):
            raise ValueError('还有未结束的分析，请先完成或在工作站取消任务：'+p.parent.name)


def prepare(source, root, windows_root, distro, port):
    source, root = Path(source).resolve(), Path(root).resolve()
    if not 1024 <= port <= 65535:
        raise ValueError('端口必须为 1024–65535')
    if any(c in windows_root for c in '\r\n\x00'):
        raise ValueError('Windows 路径无效')
    old = read_config(root)
    ensure_idle(root)
    identifier, files = inventory(source)
    for name in DIRECTORIES:
        (root/name).mkdir(parents=True, exist_ok=True)
    release = root/'程序'/identifier[:16]
    suffix = 0
    while release.exists():
        marker = release/'.source-manifest.json'
        if marker.is_file() and json.loads(marker.read_text()) == files and all(
                (release/name).is_file() and digest((release/name).read_bytes()) == hash_
                for name, hash_ in files.items()):
            break
        suffix += 1
        release = root/'程序'/(identifier[:16]+'-'+str(suffix))
    if not release.exists():
        stage = root/'程序'/('.staging-'+uuid.uuid4().hex)
        stage.mkdir()
        try:
            for name, hash_ in files.items():
                target = stage/name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source/name, target)
                if digest(target.read_bytes()) != hash_:
                    raise ValueError('复制时源码改变，请重试：'+name)
            (stage/'.source-manifest.json').write_text(json.dumps(files, sort_keys=True))
            stage.rename(release)
        finally:
            if stage.exists():
                shutil.rmtree(stage)
    dependency_key = digest((release/'requirements.txt').read_bytes())[:12]
    runtime_base = Path(os.environ.get('RNASEQ_RUNTIME_BASE', str(Path.home()/'.local/share/mouse-rnaseq-workstation')))
    runtime = runtime_base/digest(str(root).encode())[:12]/dependency_key
    now = datetime.now(timezone.utc).isoformat()
    config = dict(old, schema=1, service_id=SERVICE, root_windows=windows_root.rstrip('\\'),
                  root_wsl=str(root), program_wsl=str(release),
                  program_windows=windows_root.rstrip('\\')+'\\程序\\'+release.name,
                  runs_wsl=str(root/'分析结果'), runs_windows=windows_root.rstrip('\\')+'\\分析结果',
                  venv_wsl=str(runtime/'venv'), runtime_wsl=str(runtime), distro=distro, port=port,
                  source_digest=identifier, created_at=old.get('created_at', now), updated_at=now)
    plan = dict(config=config, previous_digest=old.get('source_digest'))
    plan_file = root/'配置'/('pending-'+uuid.uuid4().hex+'.json')
    atomic_write(plan_file, json.dumps(plan, ensure_ascii=False, indent=2).encode('utf-8'))
    return dict(config=config, plan_file=str(plan_file))


def launchers(config):
    release = Path(config['program_wsl'])
    helpers = {'station.ps1': (release/'windows/Station.ps1').read_bytes()}
    for name, action in [('启动工作站.bat', 'Start'), ('检查环境.bat', 'Check'), ('更新工作站.bat', 'Update')]:
        pause = 'pause\r\n' if action == 'Check' else 'if not "%RNASEQ_EXIT%"=="0" pause\r\n'
        helpers[name] = ('@echo off\r\nsetlocal\r\npowershell.exe -NoLogo -NoProfile -ExecutionPolicy Bypass '
                         '-File "%~dp0station.ps1" -Action '+action+'\r\n'
                         'set "RNASEQ_EXIT=%ERRORLEVEL%"\r\n'+pause+'exit /b %RNASEQ_EXIT%\r\n').encode('ascii')
    message = (f"转录组分析工作站\r\n位置：{config['root_windows']}\r\n"
               f"每日双击：启动工作站.bat\r\n网页：http://localhost:{config['port']}\r\n"
               "原始数据、参考文件、分析结果分别保存到同名目录。\r\n"
               "每个任务的图表和报告在分析结果内；下载时可手动选择导出报告目录。\r\n"
               "Python 和 Docker 计算环境由 Ubuntu/WSL 管理，未迁移已有 WSL 或 Docker 存储。\r\n"
               "Ubuntu 也可运行：bash ~/启动转录组工作站.sh\r\n")
    helpers['使用说明.txt'] = b'\xef\xbb\xbf'+message.encode('utf-8')
    return helpers


def activate(root, plan_file):
    root, plan_file = Path(root).resolve(), Path(plan_file).resolve()
    if plan_file.parent != root/'配置' or not plan_file.name.startswith('pending-'):
        raise ValueError('部署计划路径无效')
    plan = json.loads(plan_file.read_text())
    config, old = plan['config'], read_config(root)
    if old.get('source_digest') != plan['previous_digest']:
        raise ValueError('另一部署已改变当前版本，请重试，未覆盖配置')
    if Path(config['root_wsl']) != root:
        raise ValueError('部署计划与目标目录不一致')
    ensure_idle(root)
    helpers = launchers(config)
    for name, content in helpers.items():
        target = root/name
        if target.exists() and (not target.is_file() or digest(target.read_bytes()) not in
                                {digest(content), old.get('managed_files', {}).get(name)}):
            raise ValueError('启动文件已由你修改，未覆盖，请先保留或重命名：'+str(target))
    config['managed_files'] = {name: digest(content) for name, content in helpers.items()}
    for name, content in helpers.items():
        atomic_write(root/name, content)
    atomic_write(root/'配置/workstation.json', json.dumps(config, ensure_ascii=False, indent=2).encode('utf-8'))
    plan_file.unlink()
    return config


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    stage = sub.add_parser('prepare')
    stage.add_argument('source'); stage.add_argument('root'); stage.add_argument('windows_root')
    stage.add_argument('distro'); stage.add_argument('port', type=int)
    commit = sub.add_parser('activate'); commit.add_argument('root'); commit.add_argument('plan_file')
    args = parser.parse_args()
    try:
        result = prepare(args.source, args.root, args.windows_root, args.distro, args.port) if args.action == 'prepare' else activate(args.root, args.plan_file)
        print(json.dumps(result))
    except (ValueError, OSError) as e:
        parser.exit(1, str(e)+'\n')
