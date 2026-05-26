from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional

from boltz.data.types import Record


@dataclass
class Sample:
    """A sample with optional chain and interface IDs.

    Attributes
    ----------
    record : Record
        The record.
    chain_id : Optional[int]
        The chain ID.
    interface_id : Optional[int]
        The interface ID.
    """

    record_id: str
    chain_id: Optional[int] = None
    interface_id: Optional[int] = None
    weight: Optional[float] = None


class Sampler(ABC):
    """Abstract base class for samplers."""

    @abstractmethod
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
        raise NotImplementedError
