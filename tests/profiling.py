import torch
import gc


def clear_gradients(*args):
    for arg in args:
        if isinstance(arg, torch.Tensor) and arg.grad is not None:
            arg.grad = None


def clear_memory(device):
    torch._C._cuda_clearCublasWorkspaces()
    torch._dynamo.reset()
    gc.collect()
    torch.cuda.empty_cache()
    torch.cuda.reset_peak_memory_stats(device)


def peak_memory(f, *args, device):
    for _ in range(3):
        # Clean everything
        clear_memory(device)
        clear_gradients(*args)

        # Run once
        f(*args)

        # Measure peak memory
        torch.cuda.synchronize()
        memory = torch.cuda.max_memory_allocated(device)

    return memory


def current_memory(device):
    return torch.cuda.memory_allocated(device) / (1024**3)


def memory_measure(f, device, num_iters=3):
    # Clean everything
    clear_memory(device)

    # Run measurement
    print("Current memory: ", current_memory(device))
    memory = peak_memory(f, device=device)

    print("Peak memory: ", memory / (1024**3))
    return memory / (1024**3)


def memory_measure_simple(f, device, *args, **kwargs):
    # Clean everything
    clear_memory(device)
    clear_gradients(*args)

    current = current_memory(device)

    # Run once
    out = f(*args, **kwargs)

    # Measure peak memory
    torch.cuda.synchronize()
    memory = torch.cuda.max_memory_allocated(device)
    memory = memory / (1024**3)
    memory = memory - current

    return out, memory
