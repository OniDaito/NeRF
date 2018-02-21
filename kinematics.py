"""
The NeRF algorithm moved to Tensorflow for inverse kinematics

"""
import tensorflow as tf
import numpy as np
import math, itertools

# Data constants for amino acid backbones
bond_lengths = { "N_TO_A" : 1.4615, "A_TO_C" : 1.53, "C_TO_N" : 1.325 }
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

# Test Data obtained from 3NH7_1 
test_cdr_length = 11
# For now, assume we have the correct omegas for the prolines. It just makes it easy
# Torsions are all arranged here in the correct order (prev psi, prev omega, phi)

torsions = np.array([0.0, 142.95, 173.209, 
  -147.449, 138.084, -176.98,
  -110.138, 138.08, 162.29,
  -101.068, -96.169, 167.885,
  -78.779, -44.373, 175.878,
  -136.836, 164.182, -172.224,
  -63.91, 143.817, 168.896, 
  -144.503, 158.705, 175.872,
  -96.842, 103.724, -172.34,
  -85.734, -18.138, -172.979
  -150.084, 0.0, 0.0
  ])

torsions = np.array([0.0, 142.95, 
  -147.449, 138.084,
  -110.138, 138.08,
  -101.068, -96.169,
  -78.779, -44.373,
  -136.836, 164.182,
  -63.91, 143.817,
  -144.503, 158.705,
  -96.842, 103.724,
  -85.734, -18.138,
  -150.084, 0.0
  ])


torsions = np.array(list(map(math.radians, torsions)))

# bond_angles and lengths
key = "C_TO_N"
angle = bond_angles[key]
length = bond_lengths[key]
angles = []
lengths = []
angles.append(angle)
lengths.append(length)

for i in range(0, test_cdr_length * 3):
  (angle, length, key) = next_data(key)
  angles.append(angle)
  lengths.append(length)

lengths = np.array(lengths)
angles = np.array(angles)

# Initial positions
initial_positions = np.array([ [0, -1.355, 0], [0, 0, 0], [1.4466, 0.4981, 0] ])

# Target - the actual position we are after for the final Carbon Alpha
target = np.array([-0.7506869975369328, -3.2323655926888777, -6.703822569392947 ])

def normalise(x):
  ''' Normalise a vector '''
  return x / tf.sqrt(tf.reduce_sum(tf.square(x)))

def place_atoms(initial_positions, bond_angles, torsion_angles, bond_lengths, cdr_length) :
  ''' Place all of our atoms. Based on place_atom but does the entire cdr_length. '''
  positions = initial_positions
  for i in range(0, cdr_length - 1):
    idy = i * 2 
    for j in range(0, 3):
      idx = i * 3 + j
      if j == 1:
        d = place_atom_omega(positions, bond_angles[idx], bond_lengths[idx])
      else:
        d = place_atom(positions, bond_angles[idx], torsion_angles[idy], bond_lengths[idx])
        idy += 1

      positions = tf.stack([positions[1], positions[2], d], 0)

  # Return the last Carbon Alpha - should be the middle of the final 3
  return positions[1] 

def place_atom_omega(positions, bond_angle, bond_length) :
  ab = tf.subtract(positions[1], positions[0])
  bc = tf.subtract(positions[2], positions[1])
  bcn = normalise(bc)
  R = bond_length 
  
  dx = -R * tf.cos(bond_angle)
  dy = R * tf.cos(math.radians(180)) * tf.sin(bond_angle)
  dz = R * tf.sin(math.radians(180)) * tf.sin(bond_angle)

  d = tf.stack([dx,dy,dz], 0) 
  n = tf.cross(ab,bcn)
  n = normalise(n)
  nbc = tf.cross(n,bcn)
  m = tf.stack([bcn,nbc,n], 0)
  d = tf.reduce_sum(tf.multiply(tf.expand_dims(d,-1), m), axis=0)
  d = d + positions[2]
  return d

def place_atom(positions, bond_angle, torsion_angle, bond_length) :
  ''' Given the three previous atoms, the required angles and the bond
  lengths, place the next atom. Angles are in radians, lengths in angstroms.''' 
  ab = tf.subtract(positions[1], positions[0])
  bc = tf.subtract(positions[2], positions[1])
  bcn = normalise(bc)
  R = bond_length 
  
  dx = -R * tf.cos(bond_angle)
  dy = R * tf.cos(torsion_angle) * tf.sin(bond_angle)
  dz = R * tf.sin(torsion_angle) * tf.sin(bond_angle)

  d = tf.stack([dx,dy,dz], 0) 
  n = tf.cross(ab,bcn)
  n = normalise(n)
  nbc = tf.cross(n,bcn)
  m = tf.stack([bcn,nbc,n], 0)
  d = tf.reduce_sum(tf.multiply(tf.expand_dims(d,-1), m), axis=0)
  d = d + positions[2]
  return d

def basic_error(target, positions, bond_angle, torsion_angle, bond_length):
  """ Our basic error function. Reduce the difference in positions."""
  cdr_length = test_cdr_length
  d = place_atoms(positions, bond_angle, torsion_angle, bond_length, cdr_length) 
  return tf.sqrt(tf.reduce_sum(tf.pow(tf.subtract(target, d),2)))

def angle_range(x):
  """ Move our angles back into the -180 -> +180 range. """
  z = x % 360
  if z > 180:
    return z - 360
  elif z < - 180:
    return z + 360
  #return z
  return x

def to_json(results) :
  ''' Given final angles write out the json for us. '''
  import json
  frames = []

  for result in results:
    angles = list(result)
    residues = {}
    residues["data"] = [] 
    model = {}
    model["angles"] = []
    model["name"] = "3NH7_1"
    model["residues"] = ["GLU", "ARG", "TRP", "HIS", "VAL", "ARG", "GLY", "TYR", "PHE", "ASP", "HIS"]

    for res in range(0, len(model["residues"])):
      angle = {}
      if res != 0:
        angle["phi"] = angle_range(math.degrees(float(angles[res * 2 - 1])))
      else:
        angle["phi"] = 0
      if res != len(model["residues"]) - 1:
        angle["omega"] = 180.0 #angle_range(math.degrees(float(angles[res * 3 + 1])))
        angle["psi"] = angle_range(math.degrees(float(angles[res * 2])))
      else:
        angle["omega"] = 0
        angle["psi"] = 0
      
      model["angles"].append(angle)
    
    residues["data"].append(model)
    frames.append(residues)

  with open("data_angles_test.json", 'w') as f:
    f.write(json.dumps(residues))

  print(json.dumps(residues))

if __name__ == "__main__":
  # Start the tensorflow section
  with tf.Session() as sess:
    place_angle = tf.placeholder("float", None)
    place_length = tf.placeholder("float", None)
    place_torsion = tf.placeholder("float", None)
    place_position = tf.placeholder("float", None)
    place_target = tf.placeholder("float", None)

    # The variable that stands in for the real torsions. We will optimise this
    x = tf.Variable(tf.zeros([ len(torsions) ]), "torsion")

    # choose our optimiser
    error = basic_error(place_target, place_position, place_angle, x, place_length)
    train_step = tf.train.AdagradOptimizer(0.1).minimize(error)
    #train_step = tf.train.GradientDescentOptimizer(0.1).minimize(error)
    tf.global_variables_initializer().run() 
    results = []

    # Actually perform the gradient descent
    for stepnum in range(0, 100):
      sess.run([train_step], feed_dict={place_position: initial_positions, place_angle: angles, place_length: lengths, place_torsion: torsions, place_target: target})

      # Print out the actual position and torsion angles for every Carbon Alpha
      if stepnum % 1 == 0:
        pa = place_atoms(place_position, place_angle, x, place_length, test_cdr_length)
        cpos = pa.eval(feed_dict={place_position: initial_positions, place_angle: angles, place_length: lengths, place_torsion: torsions, place_target: target})
        result = sess.run(x, feed_dict={place_position: initial_positions, place_angle: angles, place_length: lengths, place_torsion: torsions, place_target: target}) 
        print(cpos) 
        results.append(result)

    # Final angles
    pa = place_atoms(place_position, place_angle, x, place_length, test_cdr_length)
    cpos = pa.eval(feed_dict={place_position: initial_positions, place_angle: angles, place_length: lengths, place_torsion: torsions, place_target: target})
    result = sess.run(x, feed_dict={place_position: initial_positions, place_angle: angles, place_length: lengths, place_torsion: torsions, place_target: target}) 
    print("Final Placing")
    print(cpos)
    results.append(result)

    print("Printing Final Angles")
    to_json(results)

  sess.close()
