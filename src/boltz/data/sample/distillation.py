from collections.abc import Iterator

from numpy.random import RandomState

from boltz.data.sample.sampler import Sample, Sampler
from boltz.data.types import Record


class DistillationSampler(Sampler):
    """A sampler for monomer distillation data."""

    def __init__(self, small_size: int = 200, small_prob: float = 0.01) -> None:
        """Initialize the sampler.

        Parameters
        ----------
        small_size : int, optional
            The maximum size to be considered small.
        small_prob : float, optional
            The probability of sampling a small item.

        """
        self._size = small_size
        self._prob = small_prob

    def sample(self, records: list[Record], random: RandomState) -> Iterator[Sample]:
        """Sample a structure from the dataset infinitely.

        Parameters
        ----------
        records : List[Record]
            The records to sample from.
        random : RandomState
            The random state for reproducibility.

        Yields
        ------
        Sample
            A data sample.

        """
        # Remove records with invalid chains
        records = [r for r in records if r.chains[0].valid]

        # Split in small and large proteins. We assume that there is only
        # one chain per record, as is the case for monomer distillation
        small = [r for r in records if r.chains[0].num_residues <= self._size]
        large = [r for r in records if r.chains[0].num_residues > self._size]

        # Sample infinitely
        while True:
            # Sample small or large
            samples = small if random.rand() < self._prob else large

            # Sample item from the list
            index = random.randint(0, len(samples))
            yield Sample(record=samples[index])
