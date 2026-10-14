"""Data preservation, activation transactions, WSL path handling and service ownership."""
from pathlib import Path
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import subprocess
import threading
import pytest
from windows.deploy import prepare, activate, REQUIRED, read_config
from windows.ubuntu_shortcut import create
from windows import station_runtime


@pytest.fixture
def station(tmp_path, monkeypatch):
    source=tmp_path/'downloaded source'
    for name in REQUIRED:
        file=source/name;file.parent.mkdir(parents=True,exist_ok=True);file.write_text('fixture '+name)
    source.joinpath('requirements.txt').write_text('flask==3.1.3\n')
    for name in ['customer-data.csv','.env','runs/job/data.csv']:
        file=source/name;file.parent.mkdir(parents=True,exist_ok=True);file.write_text('private')
    root=tmp_path/"D 盘/转录组工作站 O'Brien $(touch WRONG)"
    monkeypatch.setenv('RNASEQ_RUNTIME_BASE',str(tmp_path/'Linux runtime'))
    return source,root


def stage(source,root):
    return prepare(source,root,r"D:\转录组分析工作站",'Ubuntu-24.04',8767)


def test_prepare_then_activate_without_copying_customer_files(station):
    source,root=station
    plan=stage(source,root)
    assert not (root/'配置/workstation.json').exists()
    program=Path(plan['config']['program_wsl'])
    assert program.parent==root/'程序'
    assert not (program/'customer-data.csv').exists() and not (program/'.env').exists() and not (program/'runs').exists()
    config=activate(root,plan['plan_file'])
    assert read_config(root)['source_digest']==config['source_digest']
    assert (root/'启动工作站.bat').read_bytes().startswith(b'@echo off\r\n')
    assert config['runs_wsl']==str(root/'分析结果') and config['port']==8767


def test_repeat_install_and_updates_preserve_data_and_previous_program(station):
    source,root=station
    first=stage(source,root);activate(root,first['plan_file'])
    data=root/'原始数据/experiment.fastq';data.write_text('USER DATA')
    result=root/'分析结果/result.csv';result.write_text('USER RESULT')
    again=stage(source,root);assert again['config']['program_wsl']==first['config']['program_wsl']
    activate(root,again['plan_file'])
    source.joinpath('app.py').write_text('updated program')
    second=stage(source,root);activate(root,second['plan_file'])
    assert Path(second['config']['program_wsl'])!=Path(first['config']['program_wsl'])
    assert (Path(first['config']['program_wsl'])/'app.py').read_text()=='fixture app.py'
    assert data.read_text()=='USER DATA' and result.read_text()=='USER RESULT'


def test_modified_program_gets_separate_release_and_modified_launcher_is_preserved(station):
    source,root=station
    initial=stage(source,root);activate(root,initial['plan_file'])
    program=Path(initial['config']['program_wsl']);program.joinpath('app.py').write_text('MY EDITS')
    next_=stage(source,root);assert next_['config']['program_wsl']!=str(program)
    wrapper=root/'启动工作站.bat';wrapper.write_text('MY LAUNCHER')
    with pytest.raises(ValueError,match='未覆盖'):activate(root,next_['plan_file'])
    assert wrapper.read_text()=='MY LAUNCHER' and program.joinpath('app.py').read_text()=='MY EDITS'
    assert read_config(root)['program_wsl']==str(program)


def test_concurrent_activation_and_running_analysis_are_blocked(station):
    source,root=station
    first=stage(source,root)
    source.joinpath('app.py').write_text('next version')
    second=stage(source,root);activate(root,second['plan_file'])
    with pytest.raises(ValueError,match='另一部署'):activate(root,first['plan_file'])
    job=root/'分析结果/task';job.mkdir()
    job.joinpath('status.json').write_text(json.dumps({'state':'running'}))
    with pytest.raises(ValueError,match='未结束'):stage(source,root)


def test_ubuntu_entry_handles_unicode_shell_characters_and_tracks_current_version(station,tmp_path):
    source,root=station
    source.joinpath('windows/wsl_entry.sh').write_text('#!/usr/bin/env bash\npython3 -c \'import sys,json;print(json.dumps(sys.argv[1:]))\' "$@"\n')
    first=stage(source,root);activate(root,first['plan_file'])
    home=tmp_path/'Ubuntu home';home.mkdir()
    entry=create(root,home)
    args=json.loads(subprocess.check_output(['bash',str(entry)],text=True))
    assert args[1]==str(root) and args[2]==first['config']['program_wsl']
    assert not (tmp_path/'WRONG').exists()
    source.joinpath('app.py').write_text('updated source')
    second=stage(source,root);activate(root,second['plan_file'])
    assert create(root,home)==entry
    args=json.loads(subprocess.check_output(['bash',str(entry)],text=True))
    assert args[2]==second['config']['program_wsl']
    entry.write_text('USER SHORTCUT')
    replacement=create(root,home)
    assert entry.read_text()=='USER SHORTCUT' and replacement!=entry


@pytest.fixture
def health_server():
    payload={'service_id':'other-service'}
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200);self.end_headers();self.wfile.write(json.dumps(payload).encode())
        def log_message(self,*args):pass
    server=ThreadingHTTPServer(('127.0.0.1',0),Handler)
    thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
    yield payload,server.server_port
    server.shutdown();thread.join();server.server_close()


def test_service_identity_and_pid_reuse_do_not_stop_other_programs(tmp_path,health_server,monkeypatch):
    payload,port=health_server
    calls=[];monkeypatch.setattr(station_runtime.os,'kill',lambda *args:calls.append(args))
    assert station_runtime.probe(tmp_path,port)['state']=='occupied'
    with pytest.raises(ValueError,match='其他程序'):station_runtime.stop(tmp_path,port)
    payload.update(service_id='mouse-bulk-rnaseq',run_dir=str(tmp_path/'分析结果'),program_dir=str(tmp_path/'program'),pid=99999999)
    assert station_runtime.probe(tmp_path,port)['state']=='running'
    with pytest.raises(ValueError,match='进程与记录不一致'):station_runtime.stop(tmp_path,port)
    assert not calls


def test_active_analysis_is_never_stopped(tmp_path,health_server,monkeypatch):
    payload,port=health_server
    payload.update(service_id='mouse-bulk-rnaseq',run_dir=str(tmp_path/'分析结果'),program_dir=str(tmp_path/'program'),pid=2)
    job=tmp_path/'分析结果/task';job.mkdir(parents=True);job.joinpath('status.json').write_text('{"state":"running"}')
    calls=[];monkeypatch.setattr(station_runtime.os,'kill',lambda *args:calls.append(args))
    with pytest.raises(ValueError,match='未结束'):station_runtime.stop(tmp_path,port)
    assert not calls
