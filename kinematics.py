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

cdr_length = 4

def next_data(key):
  ''' Loop over our bond_angles and bond_lengths '''
  ff = itertools.cycle(bond_order)
  for item in ff:
    if item == key:
      next_key = next(ff)
      break
  return (bond_angles[next_key], bond_lengths[next_key], next_key)

# For now, assume we have the correct omegas for the prolines. It just makes it easy
# Torsions are all arranged here in the correct order (prev psi, prev omega, phi)
torsions = np.array([132.7336, 0.68524, -109.8387, 4.547, -161.6366, -98.671, -24.497, -171.446, -131.583])
torsions = np.array(list(map(math.radians, torsions)))

# bond_angles and lengths
key = "C_TO_N"
angle = bond_angles[key]
length = bond_lengths[key]

angles = []
lengths = []
angles.append(angle)
lengths.append(length)

for torsion in torsions:
  (angle, length, key) = next_data(key)
  angles.append(angle)
  lengths.append(length)

lengths = np.array(lengths)
angles = np.array(angles)

# Initial positions
positions = np.array([[0, -1.355, 0], [0, 0, 0], [1.4466, 0.4981, 0]]) 

# Target - the actual position we are after
target = np.array([ 1.71076135, 1.45090177, -0.88205021])

# Start the tensorflow section
with tf.Session() as sess:

  place_angle = tf.placeholder("float", None)
  place_length = tf.placeholder("float", None)
  place_torsion = tf.placeholder("float", None)
  place_position = tf.placeholder("float", None)
  place_target = tf.placeholder("float", None)

  def normalise(x):
    return x / tf.sqrt(tf.reduce_sum(tf.square(x)))
  
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
    d = tf.stack([dx,dy,dz],0) 
    n = tf.cross(ab,bcn)
    n = normalise(n)
    nbc = tf.cross(n,bcn)
    m = tf.stack([bcn,nbc,n],1)
    #https://cmsdk.com/python/how-can-i-multiply-a-vector-and-a-matrix-in-tensorflow-without-reshaping.html
    d = tf.reduce_sum(tf.multiply(tf.expand_dims(d,-1), m), axis=0)
    d = d + positions[2]
    return d

  # The variable that stands in for the real torsions. We will optimise this
  x = tf.Variable(tf.zeros([9]),"torsion")
  sess.run(x.initializer)

  def basic_error(target, positions, bond_angle, torsion_angle, bond_length):
    """ Our basic error function. Reduce the difference in positions."""
    d = place_atom(positions, bond_angle,torsion_angle, bond_length) 
    return tf.sqrt(tf.reduce_sum(tf.pow(tf.subtract(target, d),2)))

  # choose our optimiser
  error = basic_error(place_target, place_position, place_angle[0], x[0], place_length[0])
  train_step = tf.train.GradientDescentOptimizer(0.01).minimize(error)
  
  # Actually perform the gradient descent
  for stepnum in range(0,300):
    sess.run([train_step], feed_dict={place_position: positions, place_angle: angles, place_length: lengths, place_torsion: torsions, place_target: target})

    # Print out the actual position and torsion angle
    if stepnum % 10 == 0:
      pa = place_atom(place_position, place_angle[0], x[0], place_length[0])
      cpos = pa.eval(feed_dict={place_position: positions, place_angle: angles, place_length: lengths, place_torsion: torsions, place_target: target})
      angle = sess.run(x,feed_dict={place_position: positions, place_angle: angles, place_length: lengths, place_torsion: torsions, place_target: target}) 
      print(cpos, math.degrees(angle[0]))

sess.close()
