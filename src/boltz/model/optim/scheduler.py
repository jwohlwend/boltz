import torch


class AlphaFoldLRScheduler(torch.optim.lr_scheduler._LRScheduler):
    """Implements the learning rate schedule defined AF3.

    A linear warmup is followed by a plateau at the maximum
    learning rate and then exponential decay. Note that the
    initial learning rate of the optimizer in question is
    ignored; use this class' base_lr parameter to specify
    the starting point of the warmup.

    """

    def __init__(
        self,
        optimizer: torch.optim.Optimizer,
        last_epoch: int = -1,
        base_lr: float = 0.0,
        max_lr: float = 1.8e-3,
        warmup_no_steps: int = 1000,
        start_decay_after_n_steps: int = 50000,
        decay_every_n_steps: int = 50000,
        decay_factor: float = 0.95,
    ) -> None:
        """Initialize the learning rate scheduler.

        Parameters
        ----------
        optimizer : torch.optim.Optimizer
            The optimizer.
        last_epoch : int, optional
            The last epoch, by default -1
        base_lr : float, optional
            The base learning rate, by default 0.0
        max_lr : float, optional
            The maximum learning rate, by default 1.8e-3
        warmup_no_steps : int, optional
            The number of warmup steps, by default 1000
        start_decay_after_n_steps : int, optional
            The number of steps after which to start decay, by default 50000
        decay_every_n_steps : int, optional
            The number of steps after which to decay, by default 50000
        decay_factor : float, optional
            The decay factor, by default 0.95

        """
        step_counts = {
            "warmup_no_steps": warmup_no_steps,
            "start_decay_after_n_steps": start_decay_after_n_steps,
        }

        for k, v in step_counts.items():
            if v < 0:
                msg = f"{k} must be nonnegative"
                raise ValueError(msg)

        if warmup_no_steps > start_decay_after_n_steps:
            msg = "warmup_no_steps must not exceed start_decay_after_n_steps"
            raise ValueError(msg)

        self.optimizer = optimizer
        self.last_epoch = last_epoch
        self.base_lr = base_lr
        self.max_lr = max_lr
        self.warmup_no_steps = warmup_no_steps
        self.start_decay_after_n_steps = start_decay_after_n_steps
        self.decay_every_n_steps = decay_every_n_steps
        self.decay_factor = decay_factor

        super().__init__(optimizer, last_epoch=last_epoch)

    def state_dict(self) -> dict:
        state_dict = {k: v for k, v in self.__dict__.items() if k not in ["optimizer"]}
        return state_dict

    def load_state_dict(self, state_dict):
        self.__dict__.update(state_dict)

    def get_lr(self):
        if not self._get_lr_called_within_step:
            msg = (
                "To get the last learning rate computed by the scheduler, use "
                "get_last_lr()"
            )
            raise RuntimeError(msg)

        step_no = self.last_epoch

        if step_no <= self.warmup_no_steps:
            lr = self.base_lr + (step_no / self.warmup_no_steps) * self.max_lr
        elif step_no > self.start_decay_after_n_steps:
            steps_since_decay = step_no - self.start_decay_after_n_steps
            exp = (steps_since_decay // self.decay_every_n_steps) + 1
            lr = self.max_lr * (self.decay_factor**exp)
        else:  # plateau
            lr = self.max_lr

        return [lr for group in self.optimizer.param_groups]
