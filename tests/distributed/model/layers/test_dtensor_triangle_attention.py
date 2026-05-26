# SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: MIT
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.


import warnings

import pytest
import torch
from torch.distributed.tensor import Replicate, Shard, distribute_tensor

from boltz.distributed.comm import Ring2DCommTriAttn
from boltz.distributed.manager import DistributedManager
from boltz.distributed.model.layers.triangular_attention import (
    TriangleAttentionEndingNode as DistributedTriangleAttentionEndingNode,
)
from boltz.distributed.model.layers.triangular_attention import (
    TriangleAttentionStartingNode as DistributedTriangleAttentionStartingNode,
)
from boltz.distributed.model.layers.triangular_attention import (
    _Mode,
    _RingMultiHeadTriangleAttentionImpl,
    can_run_cueq_triattn_sm100f,
    cueq_is_installed,
)
from boltz.distributed.model.modules.utils import TriAttnBackend
from boltz.model.layers.triangular_attention.attention import (
    TriangleAttentionEndingNode,
    TriangleAttentionStartingNode,
)
from boltz.testing.utils import (
    assert_all_identical,
    assert_no_percentile_upshift,
    assert_tensors_identical,
    get_param_by_key,
    init_module_params_uniform,
    init_tensors_uniform,
    seed_by_rank,
    spawn_multiprocessing,
)


def parallel_assert_triangle_attention(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    dtype,
    c_in,
    c_hidden,
    no_heads,
    mode,
    layer_state_dict,
    input_x_global_host,
    mask_global_host,
    output_expected_global_host,
    d_output_expected_global_host,
    d_input_x_expected_global_host,
    grad_params_expected_global_host,
    output_global_fp32_host: torch.Tensor | None = None,
    d_input_x_global_fp32_host: torch.Tensor | None = None,
    grad_params_fp32_global_host: dict[str, torch.Tensor] | None = None,
):
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    if torch.finfo(dtype).resolution < torch.finfo(output_expected_global_host.dtype).resolution:
        raise ValueError(
            f"Target dtype {dtype} has higher precision than reference output's dtype {output_expected_global_host.dtype}"
        )

    if ((output_global_fp32_host is None) != (d_input_x_global_fp32_host is None)) or (
        (output_global_fp32_host is not None) != (grad_params_fp32_global_host is not None)
    ):
        raise ValueError(
            "output_global_fp32_host, d_input_x_global_fp32_host, and grad_params_fp32_global_host must be either all None or all not None"
        )

    check_error_hist = output_global_fp32_host is not None

    layout_map = manager.layout_subgroups["cp"]

    # Set up communication based on mode
    if mode == _Mode.Starting:
        axis_cp = 1
    elif mode == _Mode.Ending:
        axis_cp = 0
    else:
        raise ValueError(f"Invalid mode {mode}")

    ring_comm = Ring2DCommTriAttn(manager.group["cp"], layout_map, axis_cp)

    # Create reference serial module
    dtype_to_inf = {torch.float32: 1e9, torch.float64: 1e18}
    if mode == _Mode.Starting:
        module_serial = TriangleAttentionStartingNode(c_in, c_hidden, no_heads, inf=dtype_to_inf[dtype])
    elif mode == _Mode.Ending:
        module_serial = TriangleAttentionEndingNode(c_in, c_hidden, no_heads, inf=dtype_to_inf[dtype])
    else:
        raise ValueError(f"Invalid mode {mode}")

    module_serial = module_serial.to(dtype=dtype, device=manager.device)
    module_serial.load_state_dict(layer_state_dict)

    # Create distributed module
    if mode == _Mode.Starting:
        module = DistributedTriangleAttentionStartingNode(module_serial, manager.device_mesh_subgroups, ring_comm)
    elif mode == _Mode.Ending:
        module = DistributedTriangleAttentionEndingNode(module_serial, manager.device_mesh_subgroups, ring_comm)
    else:
        raise ValueError(f"Invalid mode {mode}")

    module = module.to(device=manager.device).train()

    # Input tensors have the same sharding pattern:
    # x: (B, N, N, D) - sharded on dims 1 and 2 (N and N)
    # mask: (B, N, N) - sharded on dims 1 and 2 (N and N)
    placements = (Shard(0), Shard(1), Shard(2))

    # Distribute input tensors
    input_x_dtensor = distribute_tensor(
        input_x_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    ).requires_grad_(True)

    mask_dtensor = distribute_tensor(
        mask_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )

    # Distribute expected outputs
    d_output_expected_dtensor = distribute_tensor(
        d_output_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )
    output_expected_dtensor = distribute_tensor(
        output_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )
    d_input_x_expected_dtensor = distribute_tensor(
        d_input_x_expected_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
        src_data_rank=None,
    )

    # Create copies to verify inputs aren't modified
    input_x_dtensor_copy = input_x_dtensor.detach().clone().requires_grad_(True)
    mask_dtensor_copy = mask_dtensor.detach().clone()

    if check_error_hist:
        # Forward and backward pass for error histogram checking
        output_dtensor_result = module(input_x_dtensor, mask_dtensor)
        output_dtensor_result.backward(d_output_expected_dtensor)

        output_fp32_dtensor = distribute_tensor(
            output_global_fp32_host.to(device=manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements,
            src_data_rank=None,
        )

        d_input_x_fp32_dtensor = distribute_tensor(
            d_input_x_global_fp32_host.to(device=manager.device),
            device_mesh=manager.device_mesh_subgroups,
            placements=placements,
            src_data_rank=None,
        )

        assert_no_percentile_upshift(
            output_dtensor_result.to_local(),
            output_expected_dtensor.to_local(),
            output_fp32_dtensor.to_local(),
            names_input=("output_cp_fp32", "output_serial_fp64", "output_serial_fp32"),
        )

        assert_no_percentile_upshift(
            input_x_dtensor.grad.to_local(),
            d_input_x_expected_dtensor.to_local(),
            d_input_x_fp32_dtensor.to_local(),
            names_input=("d_input_x_cp_fp32", "d_input_x_serial_fp64", "d_input_x_serial_fp32"),
        )

        # Check parameter gradients error histograms
        for name, grad_param_expected_global in grad_params_expected_global_host.items():
            grad_param_result_global = get_param_by_key(module, name).grad.full_tensor().cpu()
            assert_no_percentile_upshift(
                grad_param_result_global,
                grad_param_expected_global.to(dtype=grad_param_result_global.dtype),
                grad_params_fp32_global_host[name],
                names_input=(f"d_{name}_cp_fp32", f"d_{name}_serial_fp64", f"d_{name}_serial_fp32"),
            )
    else:
        # Forward pass
        output_dtensor_result = module(input_x_dtensor, mask_dtensor)

        # Verify inputs weren't modified
        assert_tensors_identical(
            input_x_dtensor_copy.to_local(), input_x_dtensor.to_local(), check_grad=False, check_grad_fn=False
        )
        assert_tensors_identical(mask_dtensor_copy.to_local(), mask_dtensor.to_local())

        # Test forward pass results
        assert output_dtensor_result.shape == output_expected_dtensor.shape
        assert output_dtensor_result.stride() == output_expected_dtensor.stride()
        torch.testing.assert_close(output_dtensor_result.to_local(), output_expected_dtensor.to_local())

        # Backward pass
        d_output_expected_dtensor_copy = d_output_expected_dtensor.detach().clone()
        output_dtensor_result.backward(d_output_expected_dtensor)

        # Verify upstream gradient wasn't modified
        assert_tensors_identical(d_output_expected_dtensor_copy.to_local(), d_output_expected_dtensor.to_local())

        # Test input gradients
        assert input_x_dtensor.grad.shape == d_input_x_expected_dtensor.shape
        assert input_x_dtensor.grad.stride() == d_input_x_expected_dtensor.stride()
        torch.testing.assert_close(input_x_dtensor.grad.to_local(), d_input_x_expected_dtensor.to_local())

        # Test full tensor gathering - verify distributed results match serial results
        output_global_result_host = output_dtensor_result.full_tensor().cpu()
        d_input_x_global_result_host = input_x_dtensor.grad.full_tensor().cpu()

        # Verify full tensors match expected results
        torch.testing.assert_close(output_global_result_host, output_expected_global_host.to(dtype=dtype))
        torch.testing.assert_close(d_input_x_global_result_host, d_input_x_expected_global_host.to(dtype=dtype))

        # Test parameter gradients
        grad_params_result_dtensors = {}
        for name, param in module.named_parameters():
            if param.grad is not None:
                if name not in grad_params_expected_global_host:
                    # do an extra check here to make sure the parallel computation don't result in extra gradients
                    raise ValueError(f"Parameter {name} has a resulting gradient but it is not in the reference module")
                grad_params_result_dtensors[name] = param.grad

        for name, grad_param_expected_global_host in grad_params_expected_global_host.items():
            assert name in grad_params_result_dtensors, f"Parameter {name}'s gradient is not found in result gradients"
            grad_params_result = grad_params_result_dtensors[name]
            assert grad_params_result.shape == grad_param_expected_global_host.shape
            assert grad_params_result.stride() == grad_param_expected_global_host.stride()
            grad_params_result_global = grad_params_result.full_tensor()
            torch.testing.assert_close(grad_params_result_global.cpu(), grad_param_expected_global_host.to(dtype=dtype))
            assert_all_identical(grad_params_result_global, manager.group["cp"])

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env, dtype, check_error_hist",
    (
        params_test := [
            (((1, (2, 2)), True, "cuda", "ENV"), torch.float32, False),
            (((1, (2, 2)), True, "cuda", "ENV"), torch.float64, False),
            (((2, (2, 2)), True, "cuda", "ENV"), torch.float32, False),
            (((1, (3, 3)), True, "cuda", "ENV"), torch.float32, False),
            (((1, (3, 3)), True, "cpu", "ENV"), torch.float32, False),
        ]
    ),
    indirect=["setup_env"],
    ids=[
        f"dp:{x[0][0][0]}, cp:{x[0][0][1]}, specify_method:{x[0][1]}, device_type:{x[0][2]}, method_init:{x[0][3]}, "
        f"dtype:{x[1]}, check_error_hist:{x[2]}"
        for x in params_test
    ],
)
@pytest.mark.parametrize("mode", [_Mode.Starting, _Mode.Ending])
def test_triangle_attention_parallel(setup_env, dtype, check_error_hist, mode):
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    # dtype is the dtype used by the parallel computation
    # check_error_hist determine whether to compare the error histograms between
    # (CP_in_FP32, serial_in_FP64) and (serial_in_FP32, serial_in_FP64)
    # Typically, check_error_hist will use large input dimensions to emulate
    # the real-world use cases. Same with dtype==torch.float64.

    if device_type == "cuda":
        if not torch.cuda.is_available():
            pytest.skip("skip cuda test because torch.cuda.is_available == False")
        if torch.cuda.device_count() < world_size:
            pytest.skip(f"skip cuda test because torch.cuda.device_count() != {world_size}")

    if check_error_hist:
        if grid_group_sizes["dp"] > 1:
            pytest.skip("skip error histogram check for dp > 1 to save test time")

    # For float64 and error histogram check, we use a realistic model and input size
    # with heavier computation to test the numerical stability. On the other hand,
    # a smaller model and input size incur less numerical error accumulation to allow
    # a larger range of input values to detect logical bugs inexpensively by using
    # smaller dimensions.
    test_large_model = check_error_hist or dtype == torch.float64

    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    if test_large_model:
        N = size_ring * 128  # Number of tokens
        c_in = 128  # Input dimension
        c_hidden = 32  # Hidden dimension
        no_heads = 4  # Number of attention heads
        min_val_init = -5e-2 if dtype == torch.float64 else -5e-2
        max_val_init = -min_val_init
    else:
        N = size_ring * 4  # Number of tokens
        c_in = 8  # Input dimension
        c_hidden = 16  # Hidden dimension
        no_heads = 2  # Number of attention heads
        min_val_init = -0.5
        max_val_init = 0.5

    seed = 42
    seed_by_rank(0, seed=seed)

    # compute reference results with FP64
    input_x_global_fp64 = torch.empty((B, N, N, c_in), dtype=torch.float64, requires_grad=True, device=device_type)
    mask_global_fp64 = torch.randint(0, 2, (B, N, N), dtype=torch.float64, requires_grad=False, device=device_type)
    # create pure padding chunk in the mask
    mask_global_fp64[0, N // size_ring :] = 0
    mask_global_fp64[0, :, N // size_ring :] = 0

    # Create reference serial module
    if mode == _Mode.Starting:
        reference_module = TriangleAttentionStartingNode(c_in, c_hidden, no_heads, inf=1e18)
    elif mode == _Mode.Ending:
        reference_module = TriangleAttentionEndingNode(c_in, c_hidden, no_heads, inf=1e18)
    else:
        raise ValueError(f"Invalid mode {mode}")

    # Initialize parameters to ensure reproducible behavior
    # The output activation and gradient of the layer weights typically increase by 3 to 4 orders of magnitude,
    # where the ULP would be too large and numerical error distribution becomes very wide, i.e., we would have
    # very unpredictable numerical errors. That would make the test results very noisy and not very useful to
    # detect logical bugs in the code. To avoid this, we use a smaller range for the input and layer weights.
    init_tensors_uniform([input_x_global_fp64], low=min_val_init, high=max_val_init)
    init_module_params_uniform(reference_module, low=min_val_init, high=max_val_init)

    reference_module = reference_module.to(dtype=torch.float64, device=device_type).train()
    layer_state_dict_fp64 = reference_module.state_dict()

    # Run forward pass
    output_expected_global_fp64 = reference_module(input_x_global_fp64, mask_global_fp64)
    d_output_expected_global_fp64 = torch.rand_like(output_expected_global_fp64)
    output_expected_global_fp64.backward(d_output_expected_global_fp64)

    grad_params_fp64_expected_global_host = {
        name: param.grad.detach().clone().cpu() for name, param in reference_module.named_parameters()
    }

    if check_error_hist:
        input_x_global_fp32 = input_x_global_fp64.detach().clone().to(dtype=torch.float32).requires_grad_(True)
        mask_global_fp32 = mask_global_fp64.detach().clone().to(dtype=torch.float32).requires_grad_(False)

        if mode == _Mode.Starting:
            reference_module_fp32 = TriangleAttentionStartingNode(c_in, c_hidden, no_heads)
        elif mode == _Mode.Ending:
            reference_module_fp32 = TriangleAttentionEndingNode(c_in, c_hidden, no_heads)
        else:
            raise ValueError(f"Invalid mode {mode}")

        reference_module_fp32.load_state_dict(layer_state_dict_fp64)
        reference_module_fp32 = reference_module_fp32.to(dtype=torch.float32, device=device_type).train()

        output_global_fp32 = reference_module_fp32(input_x_global_fp32, mask_global_fp32)
        d_output_expected_global_fp32 = d_output_expected_global_fp64.to(dtype=torch.float32)
        output_global_fp32.backward(d_output_expected_global_fp32)

        output_global_fp32_host = output_global_fp32.detach().clone().cpu()
        d_input_x_global_fp32_host = input_x_global_fp32.grad.detach().clone().cpu()
        grad_params_fp32_global_host = {
            name: param.grad.detach().clone().cpu() for name, param in reference_module_fp32.named_parameters()
        }
    else:
        output_global_fp32_host = None
        d_input_x_global_fp32_host = None
        grad_params_fp32_global_host = None

    # Launch parallel test across all processes
    spawn_multiprocessing(
        parallel_assert_triangle_attention,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        dtype,
        c_in,
        c_hidden,
        no_heads,
        mode,
        layer_state_dict_fp64,
        input_x_global_fp64.detach().clone().cpu(),
        mask_global_fp64.detach().clone().cpu(),
        output_expected_global_fp64.detach().clone().cpu(),
        d_output_expected_global_fp64.detach().clone().cpu(),
        input_x_global_fp64.grad.detach().clone().cpu(),
        grad_params_fp64_expected_global_host,
        output_global_fp32_host,
        d_input_x_global_fp32_host,
        grad_params_fp32_global_host,
    )


# ---------------------------------------------------------------------------
# test_cueq_triattn_sm100f_util
# ---------------------------------------------------------------------------

try:
    from cuequivariance_ops_torch.triangle_attention import _can_run_sm100f

    _cueq_ops_has_can_run_sm100f = True
except ImportError:
    _cueq_ops_has_can_run_sm100f = False


@pytest.mark.parametrize(
    "dim_token, dim_hidden, is_fwd, dtype",
    [
        (8, 8, True, torch.bfloat16),
        (8, 8, False, torch.bfloat16),
        (7, 8, True, torch.bfloat16),
        (7, 8, False, torch.bfloat16),
        (8, 7, True, torch.bfloat16),
        (8, 7, False, torch.bfloat16),
        (8, 128, True, torch.bfloat16),
        (8, 128, False, torch.bfloat16),
        (8, 129, True, torch.bfloat16),
        (8, 129, False, torch.bfloat16),
        (16, 32, True, torch.bfloat16),
        (16, 32, False, torch.bfloat16),
        (16, 32, True, torch.float16),
        (16, 32, False, torch.float16),
        (16, 32, True, torch.float32),
        (16, 32, False, torch.float32),
    ],
    ids=lambda v: str(v),
)
def test_cueq_triattn_sm100f_util(dim_token, dim_hidden, is_fwd, dtype):
    """Verify can_run_cueq_triattn_sm100f matches cuEq's private _can_run_sm100f."""
    if not torch.cuda.is_available():
        pytest.skip("CUDA not available")
    if not _cueq_ops_has_can_run_sm100f:
        pytest.skip("cuequivariance_ops_torch._can_run_sm100f not available")

    device = torch.device("cuda", 0)
    q = torch.empty(1, 1, 1, dim_token, dim_hidden, device=device, dtype=dtype)
    k = torch.empty(1, 1, 1, dim_token, dim_hidden, device=device, dtype=dtype)

    cueq_result = _can_run_sm100f(q, k, training=not is_fwd)
    # _can_run_sm100f returns (can_run: bool, device_cc: list)
    expected = cueq_result[0]
    result = can_run_cueq_triattn_sm100f(device, dtype, dim_token, dim_hidden, is_fwd)
    assert result == expected, (
        f"can_run_cueq_triattn_sm100f({device}, {dtype}, {dim_token}, {dim_hidden}, {is_fwd}) = {result}, "
        f"but _can_run_sm100f(q, k, training={not is_fwd}) = {cueq_result}"
    )


# ---------------------------------------------------------------------------
# test_triangle_attention_parallel_sm100f
# ---------------------------------------------------------------------------

SM100F_BWD_WARNING_SUBSTR = "SM100f kernel expects bias to be of the same dtype as q"


def parallel_assert_sm100f_bwd_warning(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    c_in,
    c_hidden,
    no_heads,
    mode,
    layer_state_dict,
    input_x_global_host,
    mask_global_host,
    d_output_global_host,
    expect_warning,
    mock_util_always_false,
):
    """Worker: run distributed forward+backward and check SM100f warning."""
    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    if mock_util_always_false:
        import sys

        triattn_mod = sys.modules["boltz.distributed.model.layers.triangular_attention"]
        monkeypatch.setattr(triattn_mod, "can_run_cueq_triattn_sm100f", lambda *_args, **_kw: False)

    dtype = torch.bfloat16
    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    layout_map = manager.layout_subgroups["cp"]
    axis_cp = 1 if mode == _Mode.Starting else 0
    ring_comm = Ring2DCommTriAttn(manager.group["cp"], layout_map, axis_cp)

    if mode == _Mode.Starting:
        module_serial = TriangleAttentionStartingNode(c_in, c_hidden, no_heads, inf=1e9)
    else:
        module_serial = TriangleAttentionEndingNode(c_in, c_hidden, no_heads, inf=1e9)
    module_serial = module_serial.to(dtype=dtype)
    module_serial.load_state_dict(layer_state_dict)
    module_serial = module_serial.to(device=manager.device)

    if mode == _Mode.Starting:
        module = DistributedTriangleAttentionStartingNode(module_serial, manager.device_mesh_subgroups, ring_comm)
    else:
        module = DistributedTriangleAttentionEndingNode(module_serial, manager.device_mesh_subgroups, ring_comm)
    module = module.to(device=manager.device).train()

    placements = (Shard(0), Shard(1), Shard(2))
    input_x = distribute_tensor(
        input_x_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    ).requires_grad_(True)
    mask = distribute_tensor(
        mask_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )
    d_output = distribute_tensor(
        d_output_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=placements,
    )

    output = module(input_x, mask, triattn_backend=TriAttnBackend.CUEQ)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        output.backward(d_output)

    sm100f_msgs = [w for w in caught if SM100F_BWD_WARNING_SUBSTR in str(w.message)]
    if expect_warning:
        assert sm100f_msgs, f"Rank {rank}: expected SM100f bwd warning but none was emitted"
    else:
        assert not sm100f_msgs, f"Rank {rank}: SM100f bwd warning(s) emitted ({len(sm100f_msgs)}): " + "; ".join(
            str(w.message) for w in sm100f_msgs
        )

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [
        ((2, (2, 2)), True, "cuda", "ENV"),
    ],
    indirect=["setup_env"],
    ids=["dp:2-cp:2x2-cuda"],
)
@pytest.mark.parametrize(
    "use_util_to_condition_fp32_cast",
    [True, False],
    ids=["util_active", "util_mocked_false"],
)
def test_triangle_attention_parallel_sm100f(setup_env, use_util_to_condition_fp32_cast):
    """Assert SM100f backward warning fires when util returns True, absent when mocked False."""
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    if not torch.cuda.is_available():
        pytest.skip("CUDA not available")
    if torch.cuda.device_count() < world_size:
        pytest.skip(f"Need {world_size} GPUs, have {torch.cuda.device_count()}")
    if not cueq_is_installed:
        pytest.skip("cuequivariance_torch is not installed")
    device_cc = torch.cuda.get_device_capability()
    if device_cc not in ((10, 0), (10, 3)):
        pytest.skip(f"GPU compute capability {device_cc} is not SM100/SM103")

    size_ring = grid_group_sizes["cp"][0]
    B = 2 * grid_group_sizes["dp"]
    N = size_ring * 8
    c_in = 8
    c_hidden = 8
    no_heads = 2
    mode = _Mode.Starting
    seed = 42
    seed_by_rank(0, seed=seed)

    input_x_global = torch.empty((B, N, N, c_in), dtype=torch.float64, device="cuda")
    mask_global = torch.ones((B, N, N), dtype=torch.float64, device="cuda")
    init_tensors_uniform([input_x_global], low=-0.5, high=0.5)

    reference_module = TriangleAttentionStartingNode(c_in, c_hidden, no_heads, inf=1e9)
    init_module_params_uniform(reference_module, low=-0.5, high=0.5)
    reference_module = reference_module.to(dtype=torch.float64, device="cuda")
    layer_state_dict = reference_module.state_dict()

    d_output_global = torch.rand((B, N, N, c_in), dtype=torch.float64, device="cuda")

    # True: real util pre-casts bias to q.dtype -> cuEq sees correct dtype -> no warning.
    # False: util mocked -> bias cast to fp32 -> cuEq internally detects SM100f and warns.
    mock_util_always_false = not use_util_to_condition_fp32_cast
    expect_warning = mock_util_always_false

    spawn_multiprocessing(
        parallel_assert_sm100f_bwd_warning,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        c_in,
        c_hidden,
        no_heads,
        mode,
        layer_state_dict,
        input_x_global.detach().clone().cpu(),
        mask_global.detach().clone().cpu(),
        d_output_global.detach().clone().cpu(),
        expect_warning,
        mock_util_always_false,
    )


# ---------------------------------------------------------------------------
# test_triangle_attention_bf16_autocast — dtype-preservation regression
# ---------------------------------------------------------------------------


_real_logsumexp = torch.logsumexp


def _logsumexp_fp32_promotion(*args, **kwargs):
    """Wrapper that reproduces CUDA autocast's logsumexp FP32 promotion on CPU."""
    result = _real_logsumexp(*args, **kwargs)
    if result.dtype in (torch.bfloat16, torch.float16):
        result = result.to(dtype=torch.float32)
    return result


def _unfixed_tiled_softmax_update(o_chunk, lse_m_chunk, amax_chunk, o=None, lse_m=None, amax=None):
    """tiled_softmax_attention_update WITHOUT the .to(dtype=lse_m_chunk.dtype) fix.

    When torch.logsumexp is separately monkeypatched to promote BF16 → FP32
    (as CUDA autocast does), the missing cast-back causes the FP32 cascade:
    lse_m(FP32) → delta_lse → sigmoid → o(FP32).
    """
    has_amax = amax_chunk is not None
    if o is None:
        return o_chunk, lse_m_chunk, amax_chunk

    if has_amax:
        d_lse_m = lse_m - lse_m_chunk
        amax_next = torch.maximum(amax_chunk, amax)
        delta_lse = amax_chunk - amax - d_lse_m
        o = o - torch.sigmoid(delta_lse) * (o - o_chunk)
        lse_m = lse_m_chunk + torch.logsumexp(
            torch.cat([(amax - amax_next) + d_lse_m, amax_chunk - amax_next], dim=-1),
            dim=-1,
            keepdim=True,
        )
        amax = amax_next
    else:
        d_lse_m = lse_m - lse_m_chunk
        delta_lse = -d_lse_m
        o = o - torch.sigmoid(delta_lse) * (o - o_chunk)
        lse_m = lse_m - torch.nn.functional.logsigmoid(d_lse_m)
        amax = None

    return o, lse_m, amax


def parallel_assert_bf16_autocast_dtype(
    rank,
    grid_group_sizes,
    device_type,
    backend,
    env_per_rank,
    c_in,
    c_hidden,
    no_heads,
    weight_q_global_host,
    weight_k_global_host,
    weight_v_global_host,
    input_x_global_host,
    mask_global_host,
    triangle_bias_global_host,
    d_output_global_host,
    use_bf16_logsoftmax_cast,
):
    """Worker: call _RingMultiHeadTriangleAttentionImpl.apply() and assert
    output / gradient dtypes.

    torch.logsumexp is always monkeypatched to promote half → FP32 (simulating
    CUDA autocast).  The ``use_bf16_logsoftmax_cast`` flag controls whether the
    fixed tiled_softmax_attention_update (with .to(dtype=...) cast-back) or the
    unfixed version is used.
    """
    import sys

    monkeypatch = pytest.MonkeyPatch()
    if env_per_rank is not None:
        for var_name, value in env_per_rank.items():
            if value == "<INPUT_RANK>":
                monkeypatch.setenv(var_name, f"{rank}")
                continue
            monkeypatch.setenv(var_name, value)

    monkeypatch.setattr(torch, "logsumexp", _logsumexp_fp32_promotion)

    if not use_bf16_logsoftmax_cast:
        triattn_mod = sys.modules["boltz.distributed.model.layers.triangular_attention"]
        monkeypatch.setattr(triattn_mod, "tiled_softmax_attention_update", _unfixed_tiled_softmax_update)

    dtype = torch.bfloat16
    DistributedManager.initialize(grid_group_sizes, device_type=device_type, backend=backend)
    manager = DistributedManager()

    layout_map = manager.layout_subgroups["cp"]
    axis_cp = 1
    ring_comm = Ring2DCommTriAttn(manager.group["cp"], layout_map, axis_cp)

    shard_placements = (Shard(0), Shard(1), Shard(2))
    replicate_placements = tuple(Replicate() for _ in range(manager.device_mesh_subgroups.ndim))

    q_x = distribute_tensor(
        input_x_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=shard_placements,
    ).requires_grad_(True)
    kv_x = distribute_tensor(
        input_x_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=shard_placements,
    ).requires_grad_(True)
    mask = distribute_tensor(
        mask_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=shard_placements,
    )
    triangle_bias = distribute_tensor(
        triangle_bias_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=shard_placements,
    ).requires_grad_(True)
    weight_q = distribute_tensor(
        weight_q_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=replicate_placements,
    ).requires_grad_(True)
    weight_k = distribute_tensor(
        weight_k_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=replicate_placements,
    ).requires_grad_(True)
    weight_v = distribute_tensor(
        weight_v_global_host.to(dtype=dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=replicate_placements,
    ).requires_grad_(True)

    inf_val = 1e9
    with torch.amp.autocast("cpu", dtype=dtype):
        output = _RingMultiHeadTriangleAttentionImpl.apply(
            q_x,
            kv_x,
            mask,
            triangle_bias,
            weight_q,
            weight_k,
            weight_v,
            no_heads,
            c_hidden,
            ring_comm,
            inf_val,
            TriAttnBackend.REFERENCE,
        )

    if use_bf16_logsoftmax_cast:
        assert output.dtype == dtype, f"Rank {rank}: fwd output dtype {output.dtype}, expected {dtype} (fix active)"
    else:
        assert (
            output.dtype == torch.float32
        ), f"Rank {rank}: fwd output dtype {output.dtype}, expected float32 (bug should manifest)"

    d_output = distribute_tensor(
        d_output_global_host.to(dtype=output.dtype, device=manager.device),
        device_mesh=manager.device_mesh_subgroups,
        placements=shard_placements,
    )

    if use_bf16_logsoftmax_cast:
        output.backward(d_output)
        for name, tensor in [("q_x", q_x), ("kv_x", kv_x), ("triangle_bias", triangle_bias)]:
            assert (
                tensor.grad.dtype == dtype
            ), f"Rank {rank}: bwd grad {name} dtype {tensor.grad.dtype}, expected {dtype}"
    else:
        # With the bug, the FP32 output produces FP32 do_local in backward,
        # which mixes with BF16 saved tensors (q, kT, v).  CPU matmul rejects
        # this mixed-dtype operand pair, confirming the bug propagates into
        # backward.  On CUDA the promotion would silently succeed but produce
        # FP32 gradients.
        with pytest.raises(RuntimeError, match="expected scalar type"):
            output.backward(d_output)

    DistributedManager.cleanup()
    monkeypatch.undo()


@pytest.mark.parametrize(
    "setup_env",
    [((1, (3, 3)), True, "cpu", "ENV")],
    indirect=["setup_env"],
    ids=["dp:1-cp:3x3-cpu"],
)
@pytest.mark.parametrize(
    "use_bf16_logsoftmax_cast",
    [True, False],
    ids=["fixed", "buggy"],
)
def test_triangle_attention_bf16_autocast(setup_env, use_bf16_logsoftmax_cast):
    """Regression: BF16 dtype preservation through ring attention under autocast.

    Under CUDA autocast, torch.logsumexp promotes BF16 → FP32.  Without the
    fix (.to(dtype=lse_m_chunk.dtype) in tiled_softmax_attention_update), the
    FP32 cascades from lse_m → delta_lse → sigmoid → o after ≥3 ring steps,
    making the forward output and all backward gradients FP32.

    Both cases monkeypatch torch.logsumexp to promote half → FP32 (simulating
    CUDA autocast on CPU, since CPU logsumexp preserves dtype natively).

    use_bf16_logsoftmax_cast=True : fixed tiled_softmax_attention_update (with
        .to(dtype=...) cast-back) — logsumexp still promotes but the fix casts
        back, so all outputs/grads stay BF16.
    use_bf16_logsoftmax_cast=False: unfixed tiled_softmax_attention_update —
        FP32 from logsumexp cascades into o, confirming the bug path.

    Calls _RingMultiHeadTriangleAttentionImpl.apply() directly (not through the
    module wrapper) so that downstream linears cannot mask the FP32 promotion.
    """
    grid_group_sizes, world_size, device_type, backend, _, env_per_rank = setup_env

    size_ring = grid_group_sizes["cp"][0]
    B = 2
    N = size_ring * 4
    c_in = 8
    c_hidden = 16
    no_heads = 2
    seed = 42
    seed_by_rank(0, seed=seed)

    input_x_global = torch.randn(B, N, N, c_in, dtype=torch.float64) * 0.5
    mask_global = torch.ones(B, N, N, dtype=torch.float64)
    triangle_bias_global = torch.randn(B, N, N, no_heads, dtype=torch.float64) * 0.1
    weight_q_global = torch.randn(no_heads * c_hidden, c_in, dtype=torch.float64) * 0.1
    weight_k_global = torch.randn(no_heads * c_hidden, c_in, dtype=torch.float64) * 0.1
    weight_v_global = torch.randn(no_heads * c_hidden, c_in, dtype=torch.float64) * 0.1
    d_output_global = torch.randn(B, N, N, no_heads * c_hidden, dtype=torch.float64) * 0.1

    spawn_multiprocessing(
        parallel_assert_bf16_autocast_dtype,
        world_size,
        grid_group_sizes,
        device_type,
        backend,
        env_per_rank,
        c_in,
        c_hidden,
        no_heads,
        weight_q_global,
        weight_k_global,
        weight_v_global,
        input_x_global,
        mask_global,
        triangle_bias_global,
        d_output_global,
        use_bf16_logsoftmax_cast,
    )
