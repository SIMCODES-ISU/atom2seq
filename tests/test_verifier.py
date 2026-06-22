from atom2seq.parsers import parse_xyz
from atom2seq.verifier import verify_mol

prefix = __file__.removesuffix("test_verifier.py")


def mol(name):
    return parse_xyz(prefix + f"assets/verifier_tests/{name}.xyz")


def test_valid():
    is_valid, text = verify_mol(mol("valid"))
    assert is_valid
    assert text == "This is a possible peptide."


def test_invalid():
    texts = {
        "too_small": "it only has 9 atoms, and is too small to be a peptide.",  # noqa
        "invalid_element": "Cl is not an atom found in any encoded amino acid.",  # noqa
        "missingC": "it has no C atoms in it, so it cannot be a peptide.",  # noqa
        "missingH": "it has no H atoms in it, so it cannot be a peptide.",  # noqa
        "missingN": "it has no N atoms in it, so it cannot be a peptide.",  # noqa
        "missingO": "it has no O atoms in it, so it cannot be a peptide.",  # noqa
    }  # noqa
    invalid_mols = [
        "too_small",
        "invalid_element",
        "missingC",
        "missingH",
        "missingN",
        "missingO",
    ]
    for string in invalid_mols:
        is_valid, text = verify_mol(mol(string))
        assert not is_valid
        assert text == texts[string]
