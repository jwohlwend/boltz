import torch


def test_cpu_autocast_disable_must_use_active_device_type():
    """
    Regression test for CPU autocast behavior.

    Inside a CPU autocast context, disabling autocast with device_type="cuda"
    does not actually disable CPU autocast. Disabling with device_type="cpu"
    does.
    """
    with torch.autocast("cpu", dtype=torch.bfloat16):
        x = torch.randn(4, 4)
        y = torch.randn(4, 4)

        # Wrong pattern: does not disable CPU autocast
        with torch.autocast("cuda", enabled=False):
            wrong = torch.matmul(x.float(), y.float())

        # Correct pattern: disables CPU autocast
        with torch.autocast("cpu", enabled=False):
            correct = torch.matmul(x.float(), y.float())

    assert wrong.dtype == torch.bfloat16, (
        "Using device_type='cuda' should leave CPU autocast enabled in this context"
    )
    assert correct.dtype == torch.float32, (
        "Using device_type='cpu' should disable CPU autocast in this context"
    )
