"""Regression: steering_args=None must resolve to a safe default dict."""

from boltz.model.modules.utils import DEFAULT_STEERING_ARGS, resolve_steering_args


def test_resolve_steering_args_none():
    resolved = resolve_steering_args(None)
    assert resolved["fk_steering"] is False
    assert resolved["num_particles"] == 1
    assert resolved["contact_guidance_update"] is False
    # Caller must not share the module-level dict.
    assert resolved is not DEFAULT_STEERING_ARGS


def test_resolve_steering_args_partial():
    resolved = resolve_steering_args({"fk_steering": True})
    assert resolved["fk_steering"] is True
    assert resolved["num_particles"] == DEFAULT_STEERING_ARGS["num_particles"]
    assert resolved["physical_guidance_update"] is False


def test_resolve_steering_args_full_passthrough():
    full = {
        "fk_steering": True,
        "num_particles": 3,
        "fk_lambda": 4.0,
        "fk_resampling_interval": 3,
        "physical_guidance_update": True,
        "contact_guidance_update": True,
        "num_gd_steps": 16,
    }
    resolved = resolve_steering_args(full)
    assert resolved == full
