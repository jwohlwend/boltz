"""Test suite for molecule standardization in affinity prediction.

This module tests the standardize function used in affinity prediction to ensure
it correctly handles complex multi-component molecules. The tests verify that
the RDKit implicit valence calculation fix prevents runtime errors when processing
molecules with multiple disconnected components, charged species, and complex
stereochemistry.

The fix addresses a bug where LARGEST_FRAGMENT_CHOOSER.choose() was called on
molecules created with sanitize=False without first calculating implicit valences,
causing "getNumImplicitHs() called without preceding call to calcImplicitValence()"
errors for certain complex molecular structures.
"""

import unittest

from boltz.data.parse.schema import standardize


class TestMoleculeStandardization(unittest.TestCase):
    """Test cases for molecule standardization in affinity prediction."""

    def test_standardize_simple_molecule(self):
        """Test standardization of a simple organic molecule.

        This test verifies that the standardize function works correctly
        for basic molecules that don't trigger the implicit valence bug.
        """
        smiles = "CCO"  # Ethanol

        result = standardize(smiles)

        # Should return a valid SMILES string
        self.assertIsInstance(result, str)
        self.assertGreater(len(result), 0)
        # For ethanol, the standardized form should be the same or equivalent
        self.assertEqual(result, "CCO")

    def test_standardize_amino_acid(self):
        """Test standardization of an amino acid (tyrosine).

        This test uses the same molecule from the official Boltz examples
        to ensure compatibility with known working cases.
        """
        smiles = "N[C@@H](Cc1ccc(O)cc1)C(=O)O"  # Tyrosine

        result = standardize(smiles)

        # Should return a valid SMILES string
        self.assertIsInstance(result, str)
        self.assertGreater(len(result), 0)
        # Should contain the key structural elements
        self.assertIn("N", result)
        self.assertIn("O", result)

    def test_standardize_complex_multicomponent_molecule(self):
        """Test standardization of complex multi-component molecules.

        This test specifically targets the bug where complex molecules with
        multiple disconnected components, charged species, and phosphate groups
        would cause RDKit to fail with implicit valence errors. The molecule
        tested here contains:
        - Multiple disconnected components (separated by '.')
        - Charged species ([H+])
        - Complex nucleotide-like structures with phosphate groups
        - Stereochemistry markers

        Before the fix, this would raise:
        RuntimeError: Pre-condition Violation
        getNumImplicitHs() called without preceding call to calcImplicitValence()
        """
        # Complex multi-component molecule that previously caused the error
        smiles = (
            "O=Cc1ccccc1.C#CCN.NC(=O)C1=CN(C=CC1)[C@@H]1O[C@H]"
            "(COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)"
            "n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O.[H+]"
        )

        # This should not raise an exception after the fix
        result = standardize(smiles)

        # Should return a valid SMILES string
        self.assertIsInstance(result, str)
        self.assertGreater(len(result), 0)
        # The result should be different from input due to standardization
        # (fragment selection, canonicalization, etc.)
        self.assertNotEqual(result, smiles)

    def test_standardize_molecule_with_charged_species(self):
        """Test standardization of molecules containing charged species.

        Charged species can be particularly problematic for implicit valence
        calculations, so this test ensures they are handled correctly.
        """
        smiles = "C[NH3+]"  # Methylammonium ion

        result = standardize(smiles)

        # Should return a valid SMILES string
        self.assertIsInstance(result, str)
        self.assertGreater(len(result), 0)

    def test_standardize_preserves_largest_fragment(self):
        """Test that standardization correctly selects the largest fragment.

        The LARGEST_FRAGMENT_CHOOSER should select the most significant
        component from multi-component molecules.
        """
        # Molecule with small and large fragments
        smiles = "O.CCCCCCCCCCCCCCCCCC"  # Water + long alkyl chain

        result = standardize(smiles)

        # Should return only the larger fragment (alkyl chain)
        self.assertIsInstance(result, str)
        self.assertGreater(len(result), 0)
        # Should not contain the water molecule
        self.assertNotIn("O.", result)
        self.assertNotIn(".O", result)


if __name__ == "__main__":
    unittest.main()
