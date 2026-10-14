"""Check service ownership and stop only an idle service belonging to this station."""
from pathlib import Path
import argparse
import json
import os
import signal
import socket
import time
import urllib.error
import urllib.request


def probe(root, port):
    root = Path(root).resolve()
    try:
        with urllib.request.urlopen(f'http://127.0.0.1:{port}/api/health', timeout=3) as r:
            health = json.load(r)
        owned = health.get('service_id') == 'mouse-bulk-rnaseq' and Path(health.get('run_dir', '/')).resolve() == root/'分析结果'
        return dict(state='running' if owned else 'occupied', health=health if owned else {})
    except (urllib.error.URLError, TimeoutError, OSError, ValueError):
        with socket.socket() as s:
            try:
                s.bind(('127.0.0.1', port))
            except OSError:
                return dict(state='occupied', health={})
        return dict(state='offline', health={})


def stop(root, port):
    root = Path(root).resolve()
    value = probe(root, port)
    if value['state'] == 'offline':
        return value
    if value['state'] != 'running':
        raise ValueError('端口由其他程序使用，未停止它。请改用其他端口。')
    for file in (root/'分析结果').glob('*/status.json'):
        if json.loads(file.read_text()).get('state') in ('queued', 'running'):
            raise ValueError('还有未结束的分析，请先完成或在网页中取消，未停止工作站。')
    health = value['health']
    pid = health.get('pid')
    if not isinstance(pid, int) or pid <= 1:
        raise ValueError('无法验证服务进程，未停止它。请在原终端按 Ctrl+C。')
    file = Path('/proc')/str(pid)/'cmdline'
    expected = str(Path(health['program_dir']).resolve()/'app.py').encode()
    if not file.exists() or expected not in file.read_bytes().split(b'\0'):
        raise ValueError('服务进程与记录不一致，未停止它。请在原终端按 Ctrl+C。')
    os.kill(pid, signal.SIGTERM)
    for _ in range(100):
        if probe(root, port)['state'] == 'offline':
            return dict(state='offline', health={})
        time.sleep(.1)
    raise ValueError('服务尚未退出，未强制终止，请查看原终端。')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['probe', 'wait', 'stop'])
    parser.add_argument('root'); parser.add_argument('port', type=int)
    parser.add_argument('--seconds', type=int, default=90)
    args = parser.parse_args()
    try:
        if args.action == 'wait':
            end = time.monotonic()+args.seconds
            while True:
                result = probe(args.root, args.port)
                if result['state'] == 'running':
                    break
                if result['state'] == 'occupied':
                    raise ValueError('端口已有其他程序；未打开它，请更换工作站端口。')
                if time.monotonic() >= end:
                    raise ValueError('工作站启动超时，请查看启动终端和日志/server.log。')
                time.sleep(1)
        else:
            result = stop(args.root, args.port) if args.action == 'stop' else probe(args.root, args.port)
        print(json.dumps(result))
    except (ValueError, OSError) as e:
        parser.exit(1, str(e)+'\n')
