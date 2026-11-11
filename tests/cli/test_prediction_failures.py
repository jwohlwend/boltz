"""Prediction failures must be visible to CLI callers after successful writes."""

import json
from types import SimpleNamespace
from unittest.mock import Mock

import pytest
import torch
from click.testing import CliRunner
from pytorch_lightning import LightningModule

import boltz.main as main
from boltz.data.types import Manifest
from boltz.data.write.writer import BoltzAffinityWriter, BoltzWriter
from boltz.model.models.boltz1 import Boltz1
from boltz.model.models.boltz2 import Boltz2


@pytest.fixture(autouse=True)
def _restore_torch_settings():
    enabled = torch.is_grad_enabled()
    precision = torch.get_float32_matmul_precision()
    yield
    torch.set_grad_enabled(enabled)
    torch.set_float32_matmul_precision(precision)


def _writer(tmp_path, affinity=False):
    if affinity:
        return BoltzAffinityWriter(tmp_path, tmp_path / "out")
    return BoltzWriter(tmp_path, tmp_path / "out", boltz2=True)


def _trainer(remote_failed=0):
    reduce = Mock(side_effect=lambda value, **_: value + remote_failed)
    return SimpleNamespace(strategy=SimpleNamespace(reduce=reduce))


@pytest.mark.parametrize("affinity", [False, True])
def test_failed_batch_counts_each_example_and_epoch_fails(tmp_path, affinity):
    writer = _writer(tmp_path, affinity)
    batch = {"record": [SimpleNamespace(id="one"), SimpleNamespace(id="two")]}
    writer.write_on_batch_end(None, None, {"exception": True}, [], batch, 0, 0)
    assert writer.failed == 2
    with pytest.raises(RuntimeError, match="2.*prediction.*failed"):
        writer.on_predict_epoch_end(_trainer(), SimpleNamespace(device="cpu"))


@pytest.mark.parametrize("affinity", [False, True])
def test_remote_failure_is_not_hidden_by_a_successful_local_rank(tmp_path, affinity):
    writer = _writer(tmp_path, affinity)
    trainer = _trainer(remote_failed=1)
    with pytest.raises(RuntimeError, match="1.*prediction.*failed"):
        writer.on_predict_epoch_end(trainer, SimpleNamespace(device="cpu"))
    value = trainer.strategy.reduce.call_args.args[0]
    assert value.device.type == "cpu"
    assert value.dtype == torch.int64
    assert trainer.strategy.reduce.call_args.kwargs == {"reduce_op": "sum"}


@pytest.mark.parametrize("affinity", [False, True])
def test_clean_epoch_still_participates_in_global_failure_check(tmp_path, affinity):
    writer = _writer(tmp_path, affinity)
    trainer = _trainer()
    writer.on_predict_epoch_end(trainer, SimpleNamespace(device="cpu"))
    trainer.strategy.reduce.assert_called_once()


def _run_cli(
    monkeypatch, tmp_path, structure_outcomes, affinity_outcomes=(), model="boltz2"
):
    (tmp_path / "input.yaml").write_text("version: 1\n")
    records = [
        SimpleNamespace(id=f"case{i}", affinity=bool(affinity_outcomes))
        for i in range(len(structure_outcomes))
    ]
    manifest = Manifest(records=records)
    calls = []
    monkeypatch.setattr(main, "download_boltz1", lambda *_: None)
    monkeypatch.setattr(main, "download_boltz2", lambda *_: None)
    monkeypatch.setattr(main, "check_inputs", lambda data: [data])
    monkeypatch.setattr(main, "process_inputs", lambda **_: None)
    monkeypatch.setattr(main.Manifest, "load", lambda *_: manifest)
    monkeypatch.setattr(main, "filter_inputs_structure", lambda **_: manifest)
    monkeypatch.setattr(main, "filter_inputs_affinity", lambda **_: manifest)
    monkeypatch.setattr(main, "BoltzInferenceDataModule", lambda **kw: kw)
    monkeypatch.setattr(main, "Boltz2InferenceDataModule", lambda **kw: kw)
    monkeypatch.setattr(
        main.Boltz1,
        "load_from_checkpoint",
        lambda *a, **kw: SimpleNamespace(eval=lambda: None),
    )
    monkeypatch.setattr(
        main.Boltz2,
        "load_from_checkpoint",
        lambda *a, **kw: SimpleNamespace(eval=lambda: None),
    )

    class Trainer:
        def __init__(self, callbacks, **_):
            self.callbacks = callbacks
            self.strategy = _trainer().strategy

        def predict(self, *_, **__):
            writer = self.callbacks[0]
            affinity = isinstance(writer, BoltzAffinityWriter)
            calls.append("affinity" if affinity else "structure")
            outcomes = affinity_outcomes if affinity else structure_outcomes
            for idx, failed in enumerate(outcomes):
                batch = {"record": [records[idx]]}
                if failed:
                    writer.write_on_batch_end(
                        None, None, {"exception": True}, [], batch, idx, 0
                    )
                elif affinity:
                    writer.write_on_batch_end(
                        None,
                        None,
                        {
                            "exception": False,
                            "affinity_pred_value": torch.tensor(1.25),
                            "affinity_probability_binary": torch.tensor(0.75),
                        },
                        [],
                        batch,
                        idx,
                        0,
                    )
                else:
                    dest = writer.output_dir / records[idx].id
                    dest.mkdir(parents=True, exist_ok=True)
                    (dest / "saved.cif").write_text("successful prediction\n")
            writer.on_predict_epoch_end(self, SimpleNamespace(device="cpu"))

    monkeypatch.setattr(main, "Trainer", Trainer)
    result = CliRunner().invoke(
        main.cli,
        [
            "predict",
            str(tmp_path / "input.yaml"),
            "--model",
            model,
            "--cache",
            str(tmp_path / "cache"),
            "--out_dir",
            str(tmp_path),
            "--accelerator",
            "cpu",
            "--devices",
            "1",
            "--num_workers",
            "0",
        ],
    )
    return result, calls, tmp_path / "boltz_results_input" / "predictions"


@pytest.mark.parametrize("model", ["boltz1", "boltz2"])
@pytest.mark.parametrize("outcomes", [[True], [False, True]])
def test_cli_structure_failure_is_nonzero_without_traceback(
    monkeypatch, tmp_path, model, outcomes
):
    result, calls, output = _run_cli(monkeypatch, tmp_path, outcomes, model=model)
    assert result.exit_code == 1, result.output
    assert "prediction" in result.output and "failed" in result.output
    assert isinstance(result.exception, SystemExit)
    assert calls == ["structure"]
    if outcomes[0] is False:
        assert (output / "case0" / "saved.cif").read_text() == "successful prediction\n"


def test_structure_failure_stops_affinity_before_missing_input_is_read(
    monkeypatch, tmp_path
):
    result, calls, _ = _run_cli(monkeypatch, tmp_path, [False, True], [False, False])
    assert result.exit_code == 1
    assert calls == ["structure"]
    assert "FileNotFound" not in result.output


def test_affinity_failure_keeps_successful_files_and_returns_nonzero(
    monkeypatch, tmp_path
):
    result, calls, output = _run_cli(
        monkeypatch, tmp_path, [False, False], [False, True]
    )
    assert result.exit_code == 1, result.output
    assert calls == ["structure", "affinity"]
    saved = json.loads((output / "case0" / "affinity_case0.json").read_text())
    assert saved["affinity_pred_value"] == 1.25
    assert (output / "case1" / "saved.cif").exists()
    assert not (output / "case1" / "affinity_case1.json").exists()


@pytest.mark.parametrize("affinity", [False, True])
def test_successful_cli_keeps_exit_zero(monkeypatch, tmp_path, affinity):
    result, calls, output = _run_cli(
        monkeypatch, tmp_path, [False], [False] if affinity else []
    )
    assert result.exit_code == 0, result.output
    assert calls == (["structure", "affinity"] if affinity else ["structure"])
    assert (output / "case0" / "saved.cif").exists()


@pytest.mark.parametrize("model_class", [Boltz1, Boltz2])
def test_actual_model_oom_path_reaches_failure_signal(monkeypatch, model_class):
    model = model_class.__new__(model_class)
    LightningModule.__init__(model)
    model.predict_args = dict(
        recycling_steps=1, sampling_steps=1, diffusion_samples=1, max_parallel_samples=1
    )

    def fail(*args, **kwargs):
        raise torch.cuda.OutOfMemoryError("CUDA out of memory: controlled test")

    monkeypatch.setattr(model, "forward", fail)
    assert model.predict_step({}, 0) == {"exception": True}


@pytest.mark.parametrize("model_class", [Boltz1, Boltz2])
def test_unrelated_model_error_is_not_swallowed(monkeypatch, model_class):
    model = model_class.__new__(model_class)
    LightningModule.__init__(model)
    model.predict_args = dict(
        recycling_steps=1, sampling_steps=1, diffusion_samples=1, max_parallel_samples=1
    )

    def fail(*args, **kwargs):
        raise RuntimeError("unrelated model error")

    monkeypatch.setattr(model, "forward", fail)
    with pytest.raises(RuntimeError, match="unrelated model error"):
        model.predict_step({}, 0)


@pytest.mark.parametrize("fail_rank", [-1, 1])
@pytest.mark.skipif(
    __import__("sys").platform == "win32", reason="local Gloo/fork probe"
)
def test_real_two_rank_failure_reduction(tmp_path, fail_rank):
    """Exercise the collective with a successful rank zero and failing rank one."""
    import os
    import subprocess
    import sys
    import textwrap

    script = textwrap.dedent("""
        import sys
        from pathlib import Path
        from types import SimpleNamespace
        import torch
        from torch.utils.data import DataLoader
        from pytorch_lightning import LightningModule, Trainer
        from pytorch_lightning.strategies import DDPStrategy
        from boltz.data.write.writer import BoltzAffinityWriter

        class Model(LightningModule):
            def __init__(self):
                super().__init__()
                self.weight = torch.nn.Parameter(torch.ones(1))
            def predict_step(self, batch, batch_idx, dataloader_idx=0):
                if self.global_rank == int(sys.argv[2]):
                    return {"exception": True}
                return {"exception": False, "affinity_pred_value": torch.tensor(1.25),
                        "affinity_probability_binary": torch.tensor(.75)}

        root = Path(sys.argv[1])
        writer = BoltzAffinityWriter(root, root / "out")
        trainer = Trainer(accelerator="cpu", devices=2,
            strategy=DDPStrategy(start_method="fork"), callbacks=[writer],
            logger=False, enable_checkpointing=False, enable_progress_bar=False,
            enable_model_summary=False)
        loader = DataLoader([{"record": [SimpleNamespace(id=str(i))]} for i in range(2)],
                            collate_fn=lambda batch: batch[0])
        trainer.predict(Model(), dataloaders=loader, return_predictions=False)
    """)
    env = dict(
        os.environ, CUDA_VISIBLE_DEVICES="", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1"
    )
    result = subprocess.run(
        [sys.executable, "-c", script, str(tmp_path), str(fail_rank)],
        capture_output=True,
        text=True,
        timeout=60,
        env=env,
    )
    if fail_rank < 0:
        assert result.returncode == 0, result.stdout + result.stderr
        assert (tmp_path / "out" / "1" / "affinity_1.json").exists()
    else:
        assert result.returncode != 0, result.stdout + result.stderr
        assert "1 prediction(s) failed" in result.stdout + result.stderr
        assert not (tmp_path / "out" / "1" / "affinity_1.json").exists()
    assert (
        json.loads((tmp_path / "out" / "0" / "affinity_0.json").read_text())[
            "affinity_pred_value"
        ]
        == 1.25
    )
