import time

from atom2seq.backbone_identifier import find_backbone
from atom2seq.bond_detector import bond_mol
from atom2seq.parsers import parse_cif

prefix = __file__.removesuffix("nebulin_backbone.py")

start_parse = time.time()
nebulin = parse_cif(prefix + "tests/assets/1ARK.cif")
end_parse = time.time()
print(f"Parsed 1ARK.cif in {end_parse - start_parse} seconds")

start_bond = time.time()
bond_mol(nebulin)
end_bond = time.time()
print(f"Bonded nebulin in {end_bond - start_bond} seconds")

start_bb = time.time()
find_backbone(nebulin)
end_bb = time.time()

print(f"Found nebulin's backbone in {end_bb - start_bb} seconds")
