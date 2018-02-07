"""
The NeRF algorithm

This program converts torsion angles to cartesian co-ordinates
for amino-acid back-bones. Based on the following resources:

http://onlinelibrary.wiley.com/doi/10.1002/jcc.20237/abstract
https://www.ncbi.nlm.nih.gov/pubmed/8515464
https://www.google.com/patents/WO2002073193A1?cl=en

"""

import numpy as np
import math, itertools

# TODO - PROLINE has different lengths which we should take into account
# TODO - A_TO_C angle differs by +/- 5 degrees
bond_lengths = { "N_TO_A" : 1.4615, "PRO_N_TO_A" : 1.353, "A_TO_C" : 1.53, "C_TO_N" : 1.325 }
bond_angles = { "A_TO_C" : math.radians(109), "C_TO_N" : math.radians(115), "N_TO_A" : math.radians(121) }
bond_order = ["C_TO_N", "N_TO_A", "A_TO_C"]

def next_data(key):
  ''' Loop over our bond_angles and bond_lengths '''
  ff = itertools.cycle(bond_order)
  for item in ff:
    if item == key:
      next_key = next(ff)
      break
  return (bond_angles[next_key], bond_lengths[next_key], next_key)

def place_atom(atom_a, atom_b, atom_c, bond_angle, torsion_angle, bond_length) :
  ''' Given the three previous atoms, the required angles and the bond
  lengths, place the next atom. Angles are in radians, lengths in angstroms.''' 
  # TODO - convert to sn-NeRF
  ab = np.subtract(atom_b, atom_a)
  bc = np.subtract(atom_c, atom_b)
  bcn = bc / np.linalg.norm(bc)
  R = bond_length

  # numpy is row major
  d = np.array([R * math.cos(bond_angle),
      R * math.cos(torsion_angle) * math.sin(bond_angle),
      R * math.sin(torsion_angle) * math.sin(bond_angle)])
  
  n = np.cross(ab,bcn)
  n = n / np.linalg.norm(n)
  nbc = np.cross(n,bcn)

  m = np.array([ 
    [bcn[0],nbc[0],n[0]],
    [bcn[1],nbc[1],n[1]],
    [bcn[2],nbc[2],n[2]]])

  # Still don't know why this is here :/
  d[0] = -d[0]
  d = m.dot(d)
  d = d + atom_c
  return d

if __name__ == "__main__":
  # Initial test with the first 3 atoms from 3C6S_2 
  nitrogen = np.array([0, -1.355, 0])
  alpha_carbon = np.array([0, 0, 0])
  carbon = np.array([1.4466, 0.4981, 0])          # TODO - How did we work this out again?
 
  atoms = [nitrogen, alpha_carbon, carbon]
  torsions = [132.7336, 0.68524, -109.8387, 4.547, -161.6366, -98.671, -24.497, -171.446, -131.583]
  torsions = map(math.radians, torsions)
  key = "C_TO_N"
  angle = bond_angles[key]
  length = bond_lengths[key]

  for torsion in torsions:
    atoms.append(place_atom(atoms[-3], atoms[-2], atoms[-1], angle, torsion, length))
    (angle, length, key) = next_data(key)

  print("All Atoms")
  for atom in atoms:
    print(atom)

  # Print the Carbon Alpha positions
  print("\nCarbon Alphas")
  print(atoms[1])
  print(atoms[4])
  print(atoms[7])
  print(atoms[10])

  #import atom_view
  #atom_view.launch()
