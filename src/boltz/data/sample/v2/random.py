from collections.abc import Iterator
from typing import Optional

import numpy as np
from numpy.random import RandomState

from boltz.data.sample.v2.sampler import Sampler
from boltz.data.types import Record


class RandomSampler(Sampler):
    """A simple random sampler with replacement."""

    def sample(
        self, records: list[Record], random: Optional[RandomState] = None
    ) -> Iterator[Record]:
        """Sample a structure from the dataset infinitely.

        Parameters
        ----------
        records : List[Record]
            The records to sample from.

        random : Optional[RandomState]
            Random state for reproducibility. If None, uses np.random.

        Returns
        -------
        List[Sample]
            The samples.

        """
        if random is None:
            random = np.random.default_rng()

        while True:
            item_idx = random.choice(len(records))
            yield records[item_idx]
