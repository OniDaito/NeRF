###

A short program to visualize CDR-H3 Loops

###


calculate_bond_rotation = (atom_position, previous_atom_position) ->
  difference_in_position = PXL.Math.Vec3.sub(previous_atom_position, atom_position)
  difference_in_position.normalize()
  y = new PXL.Math.Vec3(0, 1, 0)
  xp = PXL.Math.Vec3.cross(y, difference_in_position)
  dd = Math.acos(difference_in_position.dot(y))
  m = new PXL.Math.Matrix4()
  m.rotate(xp, dd)


get_test_atom_material = get_test_bond_material = () ->
  grey = new PXL.Colour.RGBA(0.8, 0.8, 0.8, 1.0)
  new PXL.Material.BasicColourMaterial(grey)


get_our_atom_material = (atom, total) ->
  progress = atom.residue.number / (total - 1)
  alpha_carbon_material_colour = new PXL.Colour.RGBA(0.8 * progress, 0 * progress, 0.4 * progress, 1.0)
  new PXL.Material.BasicColourMaterial(alpha_carbon_material_colour)


get_our_bond_material = () ->
  bond_material_colour = new PXL.Colour.RGBA(0.1, 0.8, 0.1, 1.0)
  new PXL.Material.BasicColourMaterial(bond_material_colour)


nodes_for_chain = (atoms, get_atom_material, get_bond_material) ->
  bond_geom = new PXL.Geometry.Cylinder(0.13, 50, 1, 3.82)
  atom_geom = new PXL.Geometry.Sphere(0.5, 10)

  top_node = new PXL.Node()

  for idx in [0...atoms.length]
    atom = atoms[idx]
    atom_position = atom.position
    previous_atom = atoms[idx - 1]
    # Add the atom
    atom_node = new PXL.Node(atom_geom)
    atom_node.add(get_atom_material(atom, atoms.length))
    atom_node.matrix.translate(atom_position)
    top_node.add atom_node

    if previous_atom
      # Add the bond between this atom and previous atom
      bond_node = new PXL.Node(bond_geom)
      bond_node.add get_bond_material()
      mp = PXL.Math.Vec3.add(atom_position, previous_atom.position).multScalar(0.5)
      mm = calculate_bond_rotation(atom_position, previous_atom.position)
      bond_node.matrix.translate(mp).mult(mm)
      top_node.add bond_node

  return top_node


# class ResidueOld
#   constructor : (phi, psi, omega) ->
#     @computed_alpha_carbon_position = undefined
#     # Assuming fixed bond lengths and angles with C as the central point
#     # Start left handed - initial positions
#     @a = new PXL.Math.Vec3(-2.098, 1.23, 0)  # @c_alpha
#     @b = new PXL.Math.Vec3(-1.33, 0, 0)  # @nitrogen
#     @c = new PXL.Math.Vec3(0, 0, 0)  # @carbon

#     @phi = phi
#     @psi = psi
#     @omega = omega

#     @set_positions()

#   clear_positions : () ->
#     @atom_node_a.matrix.identity()
#     @atom_node_b.matrix.identity()
#     @atom_node_b.matrix.identity()

#   set_positions : () ->
#     @clear_positions()
#     @atom_node_a.matrix.translate(@a)
#     @atom_node_b.matrix.translate(@b)
#     @atom_node_c.matrix.translate(@c)

#   # Implementation of the NeRF algorithm
#   # THIS IS THE KEY part of the program
#   # Essentially, the first 3 atoms/first residue is placed
#   # we then run next_pos which is placed, based on the previous
#   next_pos : (prev_res) ->
#     a = prev_res.a.clone()
#     b = prev_res.b.clone()
#     c = prev_res.c.clone()
#     d = @a
#     na = [@b,@c]
#     # values are taken from
#     # https://www.google.com/patents/WO2002073193A1?cl=en
#     # TODO incorporate Pro bond length of CM-NI as 1.355A
#     blengths = [1.53, 1.453, 1.325]
#     bangles = [PXL.Math.degToRad(115), PXL.Math.degToRad(109), PXL.Math.degToRad(121)]
#     torsions = [prev_res.omega, prev_res.psi, @phi]

#     for i in [0..2]
#       ab = PXL.Math.Vec3.sub(b,a)
#       abn = PXL.Math.Vec3.normalize(ab)
#       bc = PXL.Math.Vec3.sub(c,b)
#       bcn = PXL.Math.Vec3.multScalar(bc,1.0/blengths[i])
#       R = blengths[i]

#       d.x = R * Math.cos(bangles[i])
#       d.y = R * Math.cos(torsions[i]) * Math.sin(bangles[i])
#       d.z = R * Math.sin(torsions[i]) * Math.sin(bangles[i])

#       n = PXL.Math.Vec3.cross(ab,bcn).normalize()
#       nbc = PXL.Math.Vec3.cross(n,bcn)

#       #m = new PXL.Math.Matrix3([bcn.x, nbc.x, n.x, bcn.y, nbc.y, n.y, bcn.z, nbc.z, n.z])
#       m = new PXL.Math.Matrix3([bcn.x, bcn.y, bcn.z, nbc.x, nbc.y, nbc.z, n.x, n.y, n.z])
#       d.x = -d.x
#       m.multVec(d)
#       d.add(c)

#       # Shift along one
#       if i != 2
#         a = b
#         b = c
#         c = d
#         d = na[i]
#       else
#         # On the first atom which is a Ca
#         @computed_alpha_carbon_position = d

#     @set_positions()


create_chain = (model_data) ->
  alpha_carbon_positions = calculate_alpha_carbons(model_data).map (alpha_carbon) -> alpha_carbon.position
  model_node = nodes_for_chain(alpha_carbon_positions)
  return { model_node, debug: { alpha_carbon_positions } }


calculate_alpha_carbons = (model_data) ->
  backbone_atoms = calculate_backbone_atoms(model_data)
  backbone_atoms.filter (backbone_atom) -> backbone_atom.atom_type == Atom.TYPE.ALPHA_CARBON


Residue = (residue_number, data) -> ({
  phi: PXL.Math.degToRad(data.angles[residue_number].phi),
  psi: PXL.Math.degToRad(data.angles[residue_number].psi),
  omega: PXL.Math.degToRad(data.angles[residue_number].omega),
  number: residue_number,
  amino_acid: data.residues[residue_number]
})


calculate_backbone_atoms = (model_data) ->
  angles = model_data.angles
  backbone_atoms = []

  for i in [0...angles.length]
    previous_backbone_atoms = backbone_atoms[i - 1]
    residue = Residue(i, model_data)
    backbone_atoms = backbone_atoms.concat calculate_backbond_atom_positions(residue, previous_backbone_atoms, i)

  backbone_atoms

# values average taken from
# https://www.google.com/patents/WO2002073193A1?cl=en
# http://www.cryst.bbk.ac.uk/PPS95/course/3_geometry/peptide2.html
BOND_LENGTHS = {
  NITROGEN_TO_ALPHA_CARBON: (1.453 + 1.47) / 2
  PROLINE_NITROGEN_TO_ALPHA_CARBON: 1.355
  ALPHA_CARBON_TO_CARBOXYLATE_CARBON: 1.53,
  CARBOXYLATE_CARBON_TO_NITROGEN: 1.325,
}
# values average taken from
# https://www.google.com/patents/WO2002073193A1?cl=en
# TODO are Proline ALPHA_CARBON_FROM_N_TO_CC or
# bond angles different?
BOND_ANGLES = {
  VIA_ALPHA_CARBON_FROM_N_TO_CC: 109, # +/- 5
  VIA_CARBOXYLATE_CARBON_FROM_AC_TO_N: 115,
  VIA_NITROGEN_FROM_CC_TO_AC: 121
}

Atom = (x, y, z, atom_type, residue) -> ({
  position: new PXL.Math.Vec3(x, y, z),
  atom_type,
  residue
})
Atom.TYPE = {
  NITROGEN: 'NITROGEN',
  ALPHA_CARBON: 'ALPHA_CARBON',
  CARBOXYLATE_CARBON: 'CARBOXYLATE_CARBON'
}

calculate_backbond_atom_positions = (residue, previous_backbone_atom_positions, tempi) ->
  # if previous_backbone_atom_positions
  #   # do something
  # else
  nitrogen_y_position = if residue.amino_acid == 'PRO' then BOND_LENGTHS.PROLINE_NITROGEN_TO_ALPHA_CARBON else BOND_LENGTHS.NITROGEN_TO_ALPHA_CARBON
  nitrogen = Atom(0, -nitrogen_y_position, 0, Atom.TYPE.NITROGEN, residue)
  alpha_carbon = Atom(0, tempi * 2, 0, Atom.TYPE.ALPHA_CARBON, residue)

  angle = (180 - BOND_ANGLES.VIA_ALPHA_CARBON_FROM_N_TO_CC)
  hypotenuse = BOND_LENGTHS.ALPHA_CARBON_TO_CARBOXYLATE_CARBON
  carboxylate_carbon_y = Math.sin(angle) * hypotenuse
  carboxylate_carbon_x = Math.cos(angle) * hypotenuse
  carboxylate_carbon = Atom(carboxylate_carbon_x, carboxylate_carbon_y, 0, Atom.TYPE.CARBOXYLATE_CARBON, residue)

  return [
    nitrogen,
    alpha_carbon,
    carboxylate_carbon
  ]

  # Assuming fixed bond lengths and angles with C as the central point
  # Start left handed - initial positions
  # c_alpha = new PXL.Math.Vec3(-2.098, 1.23, 0)
  # nitrogen = new PXL.Math.Vec3(-1.33, 0, 0)
  # carboxylate_carbon = new PXL.Math.Vec3(0, 0, 0)

  # Implementation of the NeRF algorithm
  # Essentially, the first 3 atoms/first residue is placed
  # we then run next_pos which is placed, based on the previous
  c_alpha = prev_res.a.clone()
  nitrogen = prev_res.b.clone()
  carboxylate_carbon = prev_res.c.clone()
  d = c_alpha
  na = [nitrogen, carboxylate_carbon]

  bangles = [PXL.Math.degToRad(115), PXL.Math.degToRad(109), PXL.Math.degToRad(121)]
  torsions = [prev_res.omega, prev_res.psi, phi]

  for i in [0..2]
    ab = PXL.Math.Vec3.sub(b,a)
    abn = PXL.Math.Vec3.normalize(ab)
    bc = PXL.Math.Vec3.sub(c,b)
    bcn = PXL.Math.Vec3.multScalar(bc,1.0/blengths[i])
    R = blengths[i]

    d.x = R * Math.cos(bangles[i])
    d.y = R * Math.cos(torsions[i]) * Math.sin(bangles[i])
    d.z = R * Math.sin(torsions[i]) * Math.sin(bangles[i])

    n = PXL.Math.Vec3.cross(ab,bcn).normalize()
    nbc = PXL.Math.Vec3.cross(n,bcn)

    #m = new PXL.Math.Matrix3([bcn.x, nbc.x, n.x, bcn.y, nbc.y, n.y, bcn.z, nbc.z, n.z])
    m = new PXL.Math.Matrix3([bcn.x, bcn.y, bcn.z, nbc.x, nbc.y, nbc.z, n.x, n.y, n.z])
    d.x = -d.x
    m.multVec(d)
    d.add(c)

    # Shift along one
    if i != 2
      a = b
      b = c
      c = d
      d = na[i]
    else
      # On the first atom which is a Ca
      return d


# Actual carbon alpha positions of 3C6S_2
# from: https://www.rcsb.org/structure/3C6S
# starting at the alpha carbon of the amino acid 95 in the D chain (line 7063)
alpha_carbon_test_positions = []
alpha_carbon_test_positions.push(new PXL.Math.Vec3(31.89, 53.538, -2.462))
alpha_carbon_test_positions.push(new PXL.Math.Vec3(29.323, 54.052, -0.956))
alpha_carbon_test_positions.push(new PXL.Math.Vec3(27.71, 57.258, -2.27))
alpha_carbon_test_positions.push(new PXL.Math.Vec3(27.985, 57.642, -6.042))

# offset the test values
offset_test_positions = alpha_carbon_test_positions[0].clone()
offset_test_positions.x -= 4
test_alpha_carbons = alpha_carbon_test_positions
  .map((ac) -> PXL.Math.Vec3.sub(ac, offset_test_positions))
  .map((alpha_carbon_position) ->
    acp = alpha_carbon_position
    Atom(acp.x, acp.y, acp.z, Atom.TYPE.ALPHA_CARBON, {})
  )

# Main class for dealing with our 3D chains
class ChainsApplication

  init : () ->
    fetch("./data/data_angles_small.json")
    .then(@_parse_data)
    .then(@_setup_3d)
    .catch((err) =>
      console.error('Error initialising ', err)
    )

  draw : () ->
    # Clear and draw our shapes
    GL.clearColor(0.95, 0.95, 0.95, 1.0)
    GL.clear(GL.COLOR_BUFFER_BIT | GL.DEPTH_BUFFER_BIT)
    if @top_node
      @top_node.draw()

  _parse_data : (response) =>
    return new Promise((resolve, reject) =>
      if (response.ok)
        response.json()
        .then((data) =>
          @data = {}
          data.data.forEach((model, index) =>
            model_name_id = model.name
            if (@data[model_name_id])
              console.error('Overwriting model ' + model_name_id + ' with another model of the same name at index ' + index)
            @data[model_name_id] = model
          )
          console.log "Fetched " + Object.keys(@data).length + " models"
          resolve()
        )
        .catch((err) =>
          console.error('Error parsing data angles:', err)
          reject()
        )
      else
        console.error('Error fetching data angles: ', response.statusText)
        reject()
    )

  _error : () ->
    # Damn! Error occured
    alert("Error downloading CDR-H3 File")

  _setup_3d : () =>
    # Basic GL Functions
    GL.enable(GL.CULL_FACE)
    GL.cullFace(GL.BACK)
    GL.enable(GL.DEPTH_TEST)

    # Create the top node and add our camera
    @top_node = new PXL.Node()
    camera = new PXL.Camera.MousePerspCamera new PXL.Math.Vec3(0,0,25)
    @top_node.add camera

    # Add the test chain
    test_chain_top_node = nodes_for_chain(test_alpha_carbons, get_test_atom_material, get_test_bond_material)
    @top_node.add(test_chain_top_node)

    # For now just hard code the model to pick
    # result = create_chain(@data["3C6S_2"])
    # @top_node.add result.model_node

    uber = new PXL.GL.UberShader @top_node
    @top_node.add uber

    log_positions "Test", test_alpha_carbons
    # log_positions "Computed", result.debug.alpha_carbons


log_positions = (label, alpha_carbons) ->
    console.log(label + " carbon alpha positions")
    for a in alpha_carbons
      console.log(a.residue.amino_acid || '', ' alpha carbon at ', a.position)

chains = new ChainsApplication()

params =
  canvas : 'webgl-canvas'
  context : chains
  init : chains.init
  draw : chains.draw
  debug : true

cgl = new PXL.App params
