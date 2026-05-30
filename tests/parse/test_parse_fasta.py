import os
import pickle
from dataclasses import asdict
import pathlib
import pprint
import tempfile
import urllib.request

import pytest
import unittest

from boltz.data.parse.fasta import parse_fasta

CCD_URL = "https://huggingface.co/boltz-community/boltz-1/resolve/main/ccd.pkl"

tests_dir = os.path.dirname(os.path.abspath(__file__))
test_data_dir = os.path.join(tests_dir, 'data')

class ParseFastaTester(unittest.TestCase):

    def setUp(self):
        # Load CCD
        cache = pathlib.Path(os.path.expanduser('~/.boltz'))
        ccd = cache / "ccd.pkl"
        if not ccd.exists():
            urllib.request.urlretrieve(CCD_URL, str(ccd))  # noqa: S310

        ccd_path = ccd
        with ccd_path.open("rb") as file:
            self.ccd = pickle.load(file)  # noqa: S301

    @classmethod
    def _write_fasta(cls, lines, filename):
        with open(filename, 'w') as f:
            f.write('\n'.join(lines))

    def test_parse_fasta_onechain(self):
        fasta = [
            ">A|protein|./examples/msa/seq2.a3m",
            "QLEDSEVEAVAKGLEEMYANGVTEDNFKNYVKNNFAQQEISSVEEELNVNISDSCVANKIKDEFFAMISISAIVKAAQKKAWKELAVTVLRFAKANGLKTNAIIVAGQLALWAVQCG"
        ]
        with tempfile.NamedTemporaryFile(delete=True, suffix='.fasta') as f:
            fasta_filename = f.name
            self._write_fasta(fasta, fasta_filename)

            fasta_path = pathlib.Path(fasta_filename)
            target = parse_fasta(fasta_path, self.ccd)

            assert len(target.sequences) == 1
            sequence = target.sequences[0]
            assert len(sequence) == len(fasta[-1])
            assert str(sequence) == str(fasta[-1])
            assert target.structure.chains[0][0] == "A"

            assert target.record.chains[0].chain_name == "A"
            assert target.record.chains[0].msa_id == "./examples/msa/seq2.a3m"



    def test_parse_fasta_twochains(self):
        fasta = [
            ">A|protein",
            "MAHHHHHHVAVDAVSFTLLQDQLQSVLDTLSEREAGVVRLRFGLTDGQPRTLDEIGQVYGVTRERIRQIESKTMSKLRHPSRSQVLRDYLDGSSGSGTPEERLLRAIFGEKA",
            ">B|protein",
            "MRYAFAAEATTCNAFWRNVDMTVTALYEVPLGVCTQDPDRWTTTPDDEAKTLCRACPRRWLCARDAVESAGAEGLWAGVVIPESGRARAFALGQLRSLAERNGYPVRDHRVSAQSA"
        ]
        with tempfile.NamedTemporaryFile(delete=True, suffix='.fasta') as f:
            fasta_filename = f.name
            self._write_fasta(fasta, fasta_filename)

            fasta_path = pathlib.Path(fasta_filename)
            target = parse_fasta(fasta_path, self.ccd)

            assert len(target.sequences) == 2
            assert target.structure.chains[0][0] == "A"
            assert len(target.sequences[0]) == len(fasta[1])
            assert target.sequences[0] == fasta[1]
            assert len(target.sequences[1]) == len(fasta[3])
            assert target.sequences[1] == fasta[3]
            assert target.structure.chains[1][0] == "B"

    def test_parse_fasta_badchain(self):
        fasta = [
            ">AABB|protein",
            "QLEDSEVEAVAKGLEEMYANGVTEDNFKNYVKNNFAQQEISSVEEELNVNISDSCVANKIKDEFFAMISISAIVKAAQKKAWKELAVTVLRFAKANGLKTNAIIVAGQLALWAVQCG",
        ]
        with tempfile.NamedTemporaryFile(delete=True, suffix='.fasta') as f:
            fasta_filename = f.name
            self._write_fasta(fasta, fasta_filename)

            fasta_path = pathlib.Path(fasta_filename)

            with pytest.raises(AssertionError):
                target = parse_fasta(fasta_path, self.ccd)


