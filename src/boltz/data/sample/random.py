from collections.abc import Iterator
from dataclasses import replace

from numpy.random import RandomState

from boltz.data.sample.sampler import Sample, Sampler
from boltz.data.types import Record


class RandomSampler(Sampler):
    """A simple random sampler with replacement."""

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
        while True:
            # Sample item from the list
            index = random.randint(0, len(records))
            record = records[index]

            # Remove invalid chains and interfaces
            chains = [c for c in record.chains if c.valid]
            interfaces = [i for i in record.interfaces if i.valid]
            record = replace(record, chains=chains, interfaces=interfaces)

            yield Sample(record=record)
