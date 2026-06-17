from atom2seq.parsers import parse_xyz
from atom2seq.verifier import verify_mol

prefix = __file__.removesuffix("test_verifier.py")


def mol(name):
    return parse_xyz(prefix + f"assets/verifier_tests/{name}.xyz")


def test_valid():
    assert verify_mol(mol("valid"))


def test_invalid():
    invalid_mols = [
        "too_small",
        "invalid_element",
        "missingC",
        "missingH",
        "missingN",
        "missingO",
    ]
    for string in invalid_mols:
        assert not verify_mol(mol(string))
