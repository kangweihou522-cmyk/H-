"""Verify Docker settings preservation and rollback without touching a daemon."""
import json
import subprocess

import pytest
from scripts import configure_wsl_docker_proxy as fix

PROXY = 'http://172.20.0.1:10808'


@pytest.fixture
def engine(tmp_path, monkeypatch):
    path = tmp_path/'daemon.json'
    original = {'data-root': '/my/docker', 'proxies': {'no-proxy': 'localhost,internal.example'}, 'log-driver': 'local'}
    path.write_text(json.dumps(original))
    state = {'containers': '', 'reject_config': False, 'restarts': 0, 'fail_restart': False, 'effective': PROXY+'\n'+PROXY}
    def run(args):
        if args[-1] == '{{.ID}}': return 'local-engine'
        if args[-2:] == ['ps', '-q']: return state['containers']
        if '--validate' in args and state['reject_config']:
            raise subprocess.CalledProcessError(1, args, stderr='invalid setting')
        if args[-1] == '{{.HTTPProxy}}\n{{.HTTPSProxy}}': return state['effective']
        return ''
    def restart():
        state['restarts'] += 1
        if state['fail_restart'] and state['restarts'] == 1: raise RuntimeError('restart failed')
    monkeypatch.setattr(fix, 'run', run)
    monkeypatch.setattr(fix, 'restart_docker', restart)
    return path, original, state


def test_merge_and_backup_preserve_existing_settings(engine):
    path, original, state = engine
    fix.configure(PROXY, path)
    value = json.loads(path.read_text())
    assert value['data-root'] == original['data-root'] and value['log-driver'] == 'local'
    assert value['proxies'] == {'no-proxy': 'localhost,internal.example', 'http-proxy': PROXY, 'https-proxy': PROXY}
    backup = list(path.parent.glob('daemon.json.before-rnaseq-*'))
    assert len(backup) == 1 and json.loads(backup[0].read_text()) == original
    assert state['restarts'] == 1


@pytest.mark.parametrize('reason', ['containers', 'reject_config'])
def test_precheck_failure_does_not_change_configuration(engine, reason):
    path, original, state = engine
    state[reason] = 'running-task' if reason == 'containers' else True
    before = path.read_bytes()
    with pytest.raises((ValueError, subprocess.CalledProcessError)): fix.configure(PROXY, path)
    assert path.read_bytes() == before and state['restarts'] == 0
    assert list(path.parent.iterdir()) == [path]


@pytest.mark.parametrize('reason', ['fail_restart', 'effective'])
def test_failed_activation_restores_previous_configuration(engine, reason):
    path, original, state = engine
    state[reason] = True if reason == 'fail_restart' else ''
    before = path.read_bytes()
    with pytest.raises((RuntimeError, ValueError)): fix.configure(PROXY, path)
    assert path.read_bytes() == before and state['restarts'] == 2


def test_new_configuration_is_removed_if_activation_fails(engine):
    path, _, state = engine
    path.unlink(); state['fail_restart'] = True
    with pytest.raises(RuntimeError): fix.configure(PROXY, path)
    assert not path.exists() and state['restarts'] == 2
