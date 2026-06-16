import time

from atom2seq.backbone_identifier import find_backbone
from atom2seq.bond_detector import bond_mol
from atom2seq.parsers import parse_cif

prefix = __file__.removesuffix("titin_backbone.py")

start_parse = time.time()
titin = parse_cif(prefix + "tests/assets/1BPV.cif")
end_parse = time.time()
print(f"Parsed 1BPV.cif in {end_parse - start_parse} seconds")

start_bond = time.time()
bond_mol(titin)
end_bond = time.time()
print(f"Bonded titin in {end_bond - start_bond} seconds")

start_bb = time.time()
find_backbone(titin)
end_bb = time.time()

print(f"Found titin's backbone in {end_bb - start_bb} seconds")
