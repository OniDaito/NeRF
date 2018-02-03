###

A short program to visualize CDR-H3 Loops

###

# FILTER_ATOMS = (backbone_atom) -> backbone_atom.atom_type == Atom.TYPE.ALPHA_CARBON
FILTER_ATOMS = (backbone_atom) -> !!backbone_atom
DATA_LOCATION = "./data/data_angles.json" # "./data/data_angles_small.json"
MODEL_NAME_ID = "3C6S_2" # "4JO3_1"
RENDER_TEST_CHAIN = true


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
  progress = (atom.residue.number + 1 + total) / (2 * total)
  if atom.atom_type == Atom.TYPE.ALPHA_CARBON
    alpha_carbon_material_colour = new PXL.Colour.RGBA(0.8 * progress, 0 * progress, 0.4 * progress, 1.0)
  else if atom.atom_type == Atom.TYPE.CARBOXYLATE_CARBON
    alpha_carbon_material_colour = new PXL.Colour.RGBA(0.8 * progress, 0.8 * progress, 0.8 * progress, 1.0)
  else
    alpha_carbon_material_colour = new PXL.Colour.RGBA(0 * progress, 0.2 * progress, 0.8 * progress, 1.0)
  new PXL.Material.BasicColourMaterial(alpha_carbon_material_colour)


get_our_bond_material = () ->
  bond_material_colour = new PXL.Colour.RGBA(0.1, 0.8, 0.1, 1.0)
  new PXL.Material.BasicColourMaterial(bond_material_colour)


nodes_for_chain = (atoms, get_atom_material, get_bond_material) ->
  atom_geom = new PXL.Geometry.Sphere(0.5, 10)
  bond_geom = (length) -> new PXL.Geometry.Cylinder(0.13, 50, 1, length)

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
      bond_length = PXL.Math.Vec3.sub(atom_position, previous_atom.position).length()
      bond_position = PXL.Math.Vec3.add(atom_position, previous_atom.position).multScalar(0.5)
      bond_node = new PXL.Node(bond_geom(bond_length))
      bond_node.add get_bond_material()
      mm = calculate_bond_rotation(atom_position, previous_atom.position)
      bond_node.matrix.translate(bond_position).mult(mm)
      top_node.add bond_node

  return top_node


create_chain = (model_data) ->
  backbone_atoms = calculate_backbone_atoms(model_data).filter FILTER_ATOMS
  model_node = nodes_for_chain(backbone_atoms, get_our_atom_material, get_our_bond_material)
  return { model_node, debug: { atoms: backbone_atoms } }


Residue = (residue_number, angles, amino_acid) -> ({
  phi: PXL.Math.degToRad(angles.phi),
  psi: PXL.Math.degToRad(angles.psi),
  omega: PXL.Math.degToRad(angles.omega),
  number: residue_number,
  amino_acid: amino_acid,
  backbone_atoms: undefined # maybe mutated later
})


calculate_backbone_atoms = (model_data) ->
  residues = []
  backbone_atoms = []

  model_data.angles.forEach (angles, residue_number) ->
    residue = Residue(residue_number, angles, model_data.residues[residue_number])
    residues.push residue
    previous_residue = residues[residue_number - 1]
    backbond_atoms_for_residue = calculate_backbond_atom_positions_for_residue(residue, previous_residue)
    residue.backbone_atoms = backbond_atoms_for_residue
    backbone_atoms = backbone_atoms.concat residue.backbone_atoms.sequentially

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
PLANAR_BOND_ANGLES = {
  VIA_ALPHA_CARBON_FROM_N_TO_CC: PXL.Math.degToRad(109), # +/- 5
  VIA_CARBOXYLATE_CARBON_FROM_AC_TO_N: PXL.Math.degToRad(115),
  VIA_NITROGEN_FROM_CC_TO_AC: PXL.Math.degToRad(121)
}

Atom = ({x, y, z}, atom_type, residue) -> ({
  position: new PXL.Math.Vec3(x, y, z),
  atom_type,
  residue
})
Atom.TYPE = {
  NITROGEN: 'NITROGEN',
  ALPHA_CARBON: 'ALPHA_CARBON',
  CARBOXYLATE_CARBON: 'CARBOXYLATE_CARBON'
}


calculate_backbond_atom_positions = (last_atoms, residue, {next_atom_type, bond_length, planar_bond_angle, dihedral_bond_angle}) ->
  last_atom = last_atoms[2]
  penultimate_atom = last_atoms[1]
  antepenultimate_atom = last_atoms[0]

  a = antepenultimate_atom.position
  b = penultimate_atom.position
  c = last_atom.position

  ab = PXL.Math.Vec3.sub(b, a)
  bcn = PXL.Math.Vec3.sub(c, b).normalize()

  new_atom_position = new PXL.Math.Vec3(
    bond_length * Math.cos(planar_bond_angle),
    bond_length * Math.cos(dihedral_bond_angle) * Math.sin(planar_bond_angle),
    bond_length * Math.sin(dihedral_bond_angle) * Math.sin(planar_bond_angle),
  )

  n = PXL.Math.Vec3.cross(ab, bcn).normalize()
  nbc = PXL.Math.Vec3.cross(n, bcn)

  m = new PXL.Math.Matrix3([bcn.x, bcn.y, bcn.z, nbc.x, nbc.y, nbc.z, n.x, n.y, n.z])
  new_atom_position.x = -new_atom_position.x
  m.multVec(new_atom_position)
  new_atom_position.add(c)
  new_atom = Atom(new_atom_position, next_atom_type, residue)

  last_atoms.push new_atom
  return last_atoms.slice(1)


calculate_backbond_atom_positions_for_residue = (residue, previous_residue) ->
  nitrogen_alpha_carbon_bond_length = if residue.amino_acid == 'PRO' then BOND_LENGTHS.PROLINE_NITROGEN_TO_ALPHA_CARBON else BOND_LENGTHS.NITROGEN_TO_ALPHA_CARBON

  if previous_residue

    last_atoms = previous_residue.backbone_atoms.sequentially.slice(0)

    # Implementation of the NeRF algorithm
    [
      {
        next_atom_type: Atom.TYPE.NITROGEN,
        bond_length: BOND_LENGTHS.CARBOXYLATE_CARBON_TO_NITROGEN,
        planar_bond_angle: PLANAR_BOND_ANGLES.VIA_CARBOXYLATE_CARBON_FROM_AC_TO_N,
        dihedral_bond_angle: previous_residue.psi
      },
      {
        next_atom_type: Atom.TYPE.ALPHA_CARBON,
        bond_length: nitrogen_alpha_carbon_bond_length,
        planar_bond_angle: PLANAR_BOND_ANGLES.VIA_ALPHA_CARBON_FROM_N_TO_CC,
        dihedral_bond_angle: previous_residue.omega
      },
      {
        next_atom_type: Atom.TYPE.CARBOXYLATE_CARBON,
        bond_length: BOND_LENGTHS.ALPHA_CARBON_TO_CARBOXYLATE_CARBON,
        planar_bond_angle: PLANAR_BOND_ANGLES.VIA_CARBOXYLATE_CARBON_FROM_AC_TO_N,
        dihedral_bond_angle: residue.phi
      }
    ].forEach((params) ->
      last_atoms = calculate_backbond_atom_positions(last_atoms, residue, params)
    )

    nitrogen = last_atoms[0]
    alpha_carbon = last_atoms[1]
    carboxylate_carbon = last_atoms[2]

  else
    nitrogen = Atom({ x: 0, y: -nitrogen_alpha_carbon_bond_length, z: 0 }, Atom.TYPE.NITROGEN, residue)
    alpha_carbon = Atom({ x: 0, y: 0, z: 0 }, Atom.TYPE.ALPHA_CARBON, residue)

    angle = PXL.Math.degToRad(180) - PLANAR_BOND_ANGLES.VIA_ALPHA_CARBON_FROM_N_TO_CC
    hypotenuse = BOND_LENGTHS.ALPHA_CARBON_TO_CARBOXYLATE_CARBON
    carboxylate_carbon_position = {
      x: Math.sin(angle) * hypotenuse,
      y: Math.cos(angle) * hypotenuse,
      z: 0
    }

    carboxylate_carbon = Atom(carboxylate_carbon_position, Atom.TYPE.CARBOXYLATE_CARBON, residue)

  return {
    nitrogen,
    alpha_carbon,
    carboxylate_carbon,
    sequentially: [
      nitrogen,
      alpha_carbon,
      carboxylate_carbon
    ]
  }


get_test_alpha_carbons = (calculate_backbond_atoms) ->
  # Actual carbon alpha positions of 3C6S_2
  # from: https://www.rcsb.org/structure/3C6S
  # starting at the alpha carbon of the amino acid 95 in the D chain (line 7063)
  alpha_carbon_test_positions = []
  alpha_carbon_test_positions.push(new PXL.Math.Vec3(31.89, 53.538, -2.462))
  alpha_carbon_test_positions.push(new PXL.Math.Vec3(29.323, 54.052, -0.956))
  alpha_carbon_test_positions.push(new PXL.Math.Vec3(27.71, 57.258, -2.27))
  alpha_carbon_test_positions.push(new PXL.Math.Vec3(27.985, 57.642, -6.042))

  # offset the test values
  offset_test_position = alpha_carbon_test_positions[0].clone()
  # offset_test_position.x -= 1
  offset_test_positions = alpha_carbon_test_positions
    .map((ac) -> PXL.Math.Vec3.sub(ac, offset_test_position))

  # rotate test values round to be similar to compute values
  calculated_alpha_carbons = calculate_backbond_atoms.filter (atom) -> atom.atom_type == Atom.TYPE.ALPHA_CARBON
  rotation_vector = PXL.Math.Vec3.cross(offset_test_positions[1], calculated_alpha_carbons[1].position)
  rotation = (new PXL.Math.Matrix3())
    .rotate(rotation_vector, PXL.Math.degToRad(112))
  offset_test_positions.forEach((position) -> rotation.multVec(position) )

  # rotate aroun first vector described by two alpha carbons
  rotation = (new PXL.Math.Matrix3())
    .rotate(calculated_alpha_carbons[1].position, PXL.Math.degToRad(38))
  offset_test_positions.forEach((position) -> rotation.multVec(position) )

  test_alpha_carbons = offset_test_positions.map((alpha_carbon_position) ->
    Atom(alpha_carbon_position, Atom.TYPE.ALPHA_CARBON, {})
  )

  return test_alpha_carbons


# Main class for dealing with our 3D chains
class ChainsApplication

  init : () ->
    fetch(DATA_LOCATION)
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

    # For now just hard code the model to pick
    model_data = @data[MODEL_NAME_ID]
    if !model_data
      msg = 'No model data for ' + MODEL_NAME_ID
      alert(msg)
      throw new Error(msg)
    result = create_chain(model_data)

    # Add the test chain
    test_alpha_carbons = get_test_alpha_carbons(result.debug.atoms)
    test_chain_top_node = nodes_for_chain(test_alpha_carbons, get_test_atom_material, get_test_bond_material)

    if RENDER_TEST_CHAIN
      @top_node.add(test_chain_top_node)
    @top_node.add result.model_node

    uber = new PXL.GL.UberShader @top_node
    @top_node.add uber

    log_positions "Test", test_alpha_carbons
    log_positions "Computed", result.debug.atoms


log_positions = (label, atoms) ->
    console.log(label + " atom positions")
    for a in atoms
      console.log(a.residue.amino_acid || '', a.atom_type, 'at', a.position)

chains = new ChainsApplication()

params =
  canvas : 'webgl-canvas'
  context : chains
  init : chains.init
  draw : chains.draw
  debug : true

cgl = new PXL.App params
