#!/usr/bin/env python3
"""Isolated executable for testing launcher argument and pipe handling on Linux."""
import json
import os
import subprocess
import sys

args = sys.argv[1:]
if '--list' in args:
    sys.stderr.write('wsl: 检测到 localhost 代理配置，但未镜像到 WSL。\n')
    text = 'Ubuntu-24.04\nUbuntu\n'
    sys.stdout.buffer.write(text.encode('utf-8') if os.environ.get('WSL_UTF8') == '1' else ('\ufeff'+text).encode('utf-16-le'))
elif 'wslpath' in args:
    sys.stderr.write('wsl: 检测到 localhost 代理配置，但未镜像到 WSL。\n')
    print('/mnt/d/转录组 数据')
elif '--exec' in args:
    sys.stderr.write('wsl: 检测到 localhost 代理配置，但未镜像到 WSL。\n')
    sys.exit(subprocess.call(args[args.index('--exec') + 1:]))
elif args[0] == 'echo-argv':
    print(json.dumps(args[1:], ensure_ascii=False))
elif args[0] == 'encoding-probe':
    print(json.dumps({'WSL_UTF8': os.environ.get('WSL_UTF8')}))
    sys.stderr.write('wsl: 检测到 localhost 代理配置，但未镜像到 WSL。\n')
    sys.stderr.write('/bin/bash: line 1: --list: command not found\n')
elif args[0] == 'large-pipes':
    sys.stdout.write('o' * 200000)
    sys.stdout.flush()
    sys.stderr.write('e' * 200000)
elif args[0] == 'fail':
    sys.stderr.write('/usr/bin/w: /usr/bin/w: cannot execute binary file\n')
    sys.exit(126)
elif args[0] == 'interactive':
    # Confirm inherited stdin was not closed by the noncaptured invocation.
    if sys.stdin.read(3) != 'ABC':
        sys.exit(3)
else:
    sys.exit(2)
