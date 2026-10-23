"""Configure the local WSL Docker Engine after checking its Windows HTTP proxy."""
import argparse
import ipaddress
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import tempfile
import urllib.error
import urllib.request
import uuid

DOCKER = ['docker', '--host=unix:///var/run/docker.sock']
BASE_IMAGE = 'debian:bookworm-slim@sha256:7c7b2c966bc9ee8cedfeef67e0e279108992c77681fa595db4a9d65c06ccc587'


def run(args):
    return subprocess.run(args, check=True, text=True, capture_output=True).stdout.strip()


def windows_proxy(port):
    routes = json.loads(run(['ip', '-j', '-4', 'route', 'show', 'default']))
    gateways = [r['gateway'] for r in routes if r.get('gateway')]
    if not gateways:
        raise ValueError('未找到WSL默认网关，尚未修改Docker配置。')
    gateway = ipaddress.IPv4Address(gateways[0])
    return f'http://{gateway}:{port}'


def check_registry(proxy):
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({'http': proxy, 'https': proxy}))
    try:
        with opener.open('https://registry-1.docker.io/v2/', timeout=20) as response:
            if response.status != 200:
                raise ValueError('Docker Hub返回了意外状态。')
    except urllib.error.HTTPError as error:
        # A registry authentication challenge confirms the HTTPS connection.
        if error.code != 401 or 'auth.docker.io' not in error.headers.get('WWW-Authenticate', ''):
            raise ValueError(f'代理访问Docker Hub失败，HTTP {error.code}。') from error
    except (urllib.error.URLError, TimeoutError, OSError) as error:
        raise ValueError(f'WSL无法通过 {proxy} 访问Docker Hub：{error}。请检查代理软件允许局域网连接、端口10808及Windows防火墙；尚未修改Docker配置。') from error


def restart_docker():
    command = ['systemctl', 'restart', 'docker'] if Path('/run/systemd/system').is_dir() else ['service', 'docker', 'restart']
    run(command)
    run(DOCKER+['info', '--format', '{{.ServerVersion}}'])


def configure(proxy, config_file=Path('/etc/docker/daemon.json')):
    config_file = Path(config_file)
    if config_file.is_symlink():
        raise ValueError('daemon.json是符号链接，请先确认Docker实际配置位置；未修改。')
    if run(['docker', 'info', '--format', '{{.ID}}']) != run(DOCKER+['info', '--format', '{{.ID}}']):
        raise ValueError('当前Docker上下文与本机引擎不一致，未修改。')
    if run(DOCKER+['ps', '-q']):
        raise ValueError('Docker还有正在运行的容器。请先正常结束相关任务，再运行修复；未停止任何容器。')
    original = config_file.read_bytes() if config_file.exists() else None
    mode = stat.S_IMODE(config_file.stat().st_mode) if original is not None else 0o600
    config = json.loads(original.decode('utf-8-sig')) if original is not None else {}
    if not isinstance(config, dict) or not isinstance(config.get('proxies', {}), dict):
        raise ValueError('现有Docker配置格式异常，未修改。')
    proxies = config.setdefault('proxies', {})
    proxies.update({'http-proxy': proxy, 'https-proxy': proxy})
    proxies.setdefault('no-proxy', 'localhost,127.0.0.1,::1')
    config_file.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=config_file.parent, delete=False) as output:
        pending = Path(output.name)
        output.write(json.dumps(config, ensure_ascii=False, indent=2).encode('utf-8'))
    try:
        pending.chmod(mode)
        run(['dockerd', '--validate', '--config-file', str(pending)])
        backup = None
        if original is not None:
            backup = config_file.with_name('daemon.json.before-rnaseq-'+uuid.uuid4().hex)
            shutil.copy2(config_file, backup)
            print(f'原Docker配置备份：{backup}', flush=True)
        pending.replace(config_file)
        try:
            restart_docker()
            effective = run(DOCKER+['info', '--format', '{{.HTTPProxy}}\n{{.HTTPSProxy}}']).splitlines()
            if effective != [proxy, proxy]:
                raise ValueError('Docker没有加载新的代理配置，可能使用了自定义配置文件。')
        except Exception:
            if original is None:
                config_file.unlink(missing_ok=True)
            else:
                shutil.copy2(backup, config_file)
            print('新配置未能生效，已恢复之前的配置，正在恢复Docker服务。', flush=True)
            restart_docker()
            raise
    finally:
        pending.unlink(missing_ok=True)
    print(f'Docker代理已生效：{proxy}', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--port', type=int, default=10808)
    args = parser.parse_args()
    try:
        if not 1 <= args.port <= 65535:
            raise ValueError('代理端口无效。')
        if os.geteuid() != 0:
            raise ValueError('请通过sudo python3运行此脚本。')
        proxy = windows_proxy(args.port)
        print(f'正在检查Windows代理：{proxy}', flush=True)
        check_registry(proxy)
        print('Docker Hub连接检查通过。', flush=True)
        configure(proxy)
        print('正在下载工作站需要的基础镜像……', flush=True)
        subprocess.run(DOCKER+['pull', BASE_IMAGE], check=True)
        print('基础镜像已就绪，可以重新运行部署到D盘.bat。', flush=True)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        detail = error.stderr if isinstance(error, subprocess.CalledProcessError) else str(error)
        parser.exit(1, '配置未完成：'+str(detail or error)+'\n')
