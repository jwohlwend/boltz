from boltz.data.sample.v2.sampler import Sample, Sampler
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

    def sample(self, records: list[Record]) -> list[Sample]:
        """Sample a structure from the dataset infinitely.

        Parameters
        ----------
        records : List[Record]
            The records to sample from.

        Returns
        -------
        List[Sample]
            The samples.

        """
        # Remove records with invalid chains
        records = [r for r in records if r.chains[0].valid]

        # Split in small and large proteins. We assume that there is only
        # one chain per record, as is the case for monomer distillation
        small = [r for r in records if r.chains[0].num_residues <= self._size]
        large = [r for r in records if r.chains[0].num_residues > self._size]

        # Assign uniform weights to the proteins, with prob amount of small
        weights = [self._prob / len(small)] * len(small)
        weights += [(1 - self._prob) / len(large)] * len(large)

        # Create samples
        samples = [
            Sample(record_id=r.id, chain_id=0, weight=w)
            for r, w in zip(records, weights)
        ]
        return samples
