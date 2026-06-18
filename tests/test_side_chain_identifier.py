import pytest

from atom2seq.backbone_identifier import label_backbone
from atom2seq.bond_detector import bond_mol
from atom2seq.parsers import parse_xyz
from atom2seq.side_chain_identifier import id_side_chains

prefix = __file__.removesuffix("test_side_chain_identifier.py")


@pytest.fixture
def glycine():
    gly = parse_xyz(prefix + "assets/side_chain_identifier_tests/glycine.xyz")
    bond_mol(gly)
    label_backbone(gly)
    return gly


def test_glycine(glycine):
    assert id_side_chains(glycine) == ["G"]
