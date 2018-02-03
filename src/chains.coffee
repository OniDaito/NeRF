### ABOUT
             .__
_________  __|  |
\____ \  \/  /  |
|  |_> >    <|  |__
|   __/__/\_ \____/
|__|        \/     js
                    PXL.js
                    Benjamin Blundell - ben@pxljs.com
                    http://pxljs.com
This software is released under the MIT Licence. See LICENCE.txt for details

A short program to visualize CDR-H3 Loops

###

# Actual carbon alpha positions of 3C6S_2
# from: https://www.rcsb.org/structure/3C6S
# starting at the alpha carbon of the amino acid 95 in the D chain (line 7063)
real_ca = []
real_ca.push(new PXL.Math.Vec3(31.89, 53.538, -2.462))
real_ca.push(new PXL.Math.Vec3(29.323, 54.052, -0.956))
real_ca.push(new PXL.Math.Vec3(27.71, 57.258, -2.27))
real_ca.push(new PXL.Math.Vec3(27.985, 57.642, -6.042))

# Starting from zero, here are the real values
alpha_carbon_test_positions = []
alpha_carbon_test_positions.push(new PXL.Math.Vec3(0, 0, 0))
alpha_carbon_test_positions.push(PXL.Math.Vec3.sub(real_ca[1], real_ca[0]))
alpha_carbon_test_positions.push(PXL.Math.Vec3.sub(real_ca[2], real_ca[0]))
alpha_carbon_test_positions.push(PXL.Math.Vec3.sub(real_ca[3], real_ca[0]))


calculate_bond_rotation = (atom_position, previous_atom_position) ->
  difference_in_position = PXL.Math.Vec3.sub(previous_atom_position, atom_position)
  difference_in_position.normalize()
  y = new PXL.Math.Vec3(0, 1, 0)
  xp = PXL.Math.Vec3.cross(y, difference_in_position)
  dd = Math.acos(difference_in_position.dot(y))
  m = new PXL.Math.Matrix4()
  m.rotate(xp, dd)

nodes_for_chain = (atom_positions) ->
  bond_geom = new PXL.Geometry.Cylinder(0.13, 50, 1, 3.82)
  atom_geom = new PXL.Geometry.Sphere(0.5, 10)

  pg = new PXL.Colour.RGBA(0.8, 0.8, 0.8, 1.0)
  bond_mat = new PXL.Material.BasicColourMaterial(pg)
  atom_mat = new PXL.Material.BasicColourMaterial(pg)

  top_node = new PXL.Node()

  residue_atoms_node = new PXL.Node()
  residue_bonds_node = new PXL.Node()

  residue_atoms_node.add(atom_mat)
  residue_bonds_node.add(bond_mat)

  top_node.add residue_atoms_node
  top_node.add residue_bonds_node

  for idx in [0...atom_positions.length]
    atom_position = atom_positions[idx]
    previous_atom_position = atom_positions[idx - 1]
    # Add the atom
    atom_node = new PXL.Node(atom_geom)
    residue_atoms_node.add(atom_node)
    atom_node.matrix.translate(atom_position)

    if previous_atom_position
      # Add the bond between this atom and previous atom
      bond_node = new PXL.Node(bond_geom)
      residue_bonds_node.add(bond_node)
      mp = PXL.Math.Vec3.add(atom_position, previous_atom_position).multScalar(0.5)
      mm = calculate_bond_rotation(atom_position, previous_atom_position)
      bond_node.matrix.translate(mp).mult(mm)

  return top_node


class Residue
  constructor : (phi, psi, omega, bond_geom, atom_geom, bond_mat, atom_mat, show_bond) ->
    @computed_carbon_alpha_position = undefined
    # Assuming fixed bond lengths and angles with C as the central point
    # Start left handed - initial positions
    @a = new PXL.Math.Vec3(-2.098, 1.23, 0)  # @c_alpha
    @b = new PXL.Math.Vec3(-1.33, 0, 0)  # @nitrogen
    @c = new PXL.Math.Vec3(0, 0, 0)  # @carbon

    @phi = phi
    @psi = psi
    @omega = omega

    @residue_node = new PXL.Node()
    residue_atoms_node = new PXL.Node()
    residue_bonds_node = new PXL.Node()

    residue_atoms_node.add(atom_mat)
    residue_bonds_node.add(bond_mat)

    @atom_node_a = new PXL.Node(atom_geom)
    @atom_node_b = new PXL.Node(atom_geom)
    @atom_node_c = new PXL.Node(atom_geom)

    @bond_node_a = new PXL.Node(bond_geom)

    residue_atoms_node.add @atom_node_a
    #residue_atoms_node.add @atom_node_b
    #residue_atoms_node.add @atom_node_c

    if show_bond
      residue_bonds_node.add @bond_node_a

    @set_positions()
    # ignore bonds for now
    @residue_node.add(residue_atoms_node)
    @residue_node.add(residue_bonds_node)

  clear_positions : () ->
    @atom_node_a.matrix.identity()
    @atom_node_b.matrix.identity()
    @atom_node_b.matrix.identity()

  set_positions : () ->
    @clear_positions()
    @atom_node_a.matrix.translate(@a)
    @atom_node_b.matrix.translate(@b)
    @atom_node_c.matrix.translate(@c)

  # Implementation of the NeRF algorithm
  # THIS IS THE KEY part of the program
  # Essentially, the first 3 atoms/first residue is placed
  # we then run next_pos which is placed, based on the previous
  next_pos : (prev_res) ->
    a = prev_res.a.clone()
    b = prev_res.b.clone()
    c = prev_res.c.clone()
    d = @a
    na = [@b,@c]
    # values are taken from
    # https://www.google.com/patents/WO2002073193A1?cl=en
    # TODO incorporate Pro bond length of CM-NI as 1.355A
    blengths = [1.53, 1.453, 1.325]
    bangles = [PXL.Math.degToRad(115), PXL.Math.degToRad(109), PXL.Math.degToRad(121)]
    torsions = [prev_res.omega, prev_res.psi, @phi]

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
        @computed_carbon_alpha_position = d

    @set_positions()


get_atom_material = (i, num_residues) ->
  pink = new PXL.Colour.RGBA(0.8, 0, 0.4, 1.0)
  pg = pink.clone()
  pg.r = pg.r / num_residues * (i + 1)
  pg.g = pg.g / num_residues * (i + 1)
  pg.b = pg.b / num_residues * (i + 1)
  calpha_material = new PXL.Material.BasicColourMaterial(pg)


get_bond_material = (i, num_residues) ->
  green = new PXL.Colour.RGBA(0.1, 0.8, 0.1, 1.0)
  tg = green.clone()
  tg.r = tg.r / num_residues * (i + 1)
  tg.g = tg.g / num_residues * (i + 1)
  tg.b = tg.b / num_residues * (i + 1)
  backbone_material = new PXL.Material.BasicColourMaterial(tg)


create_chain = (model_data) ->
  residues = []
  computed_carbon_alpha_positions = []
  bond_geom = new PXL.Geometry.Cylinder(0.13, 50, 1, 3.82)
  atom_geom = new PXL.Geometry.Sphere(0.5, 10)

  model_node = new PXL.Node()
  num_residues = model_data.residues.length

  prev_res = null
  show_bond = false

  for i in [0...num_residues]
  #for i in [0..2]
    model_angles = model_data.angles[i]

    phi = PXL.Math.degToRad(model_angles.phi)
    psi = PXL.Math.degToRad(model_angles.psi)
    omega = PXL.Math.degToRad(model_angles.omega)

    residue = new Residue(phi, psi, omega, bond_geom, atom_geom, get_bond_material(i, num_residues), get_atom_material(i, num_residues), show_bond)

    if i > 0
      residue.next_pos(prev_res)

      # Now work on the bonds
      mp = PXL.Math.Vec3.add(residue.a, prev_res.a).multScalar(0.5)
      mm = calculate_bond_rotation(residue.a, prev_res.a)
      residue.bond_node_a.matrix.translate(mp).mult(mm)

    rn = residue.residue_node
    model_node.add rn
    residues.push residue
    computed_carbon_alpha_positions.push residue.computed_carbon_alpha_position

    show_bond = true
    prev_res = residue

  { model_node, debug: { computed_carbon_alpha_positions, residues } }

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
    test_chain_top_node = nodes_for_chain(alpha_carbon_test_positions)
    @top_node.add(test_chain_top_node)

    # For now just hard code the model to pick
    result = create_chain(@data["3C6S_2"])
    @top_node.add result.model_node

    uber = new PXL.GL.UberShader @top_node
    @top_node.add uber

    log_positions "Test", alpha_carbon_test_positions
    log_positions "Computed", result.debug.computed_carbon_alpha_positions


log_positions = (label, alpha_carbon_positions) ->
    console.log(label + " carbon alpha positions")
    for a in alpha_carbon_positions
      console.log(a)

chains = new ChainsApplication()

params =
  canvas : 'webgl-canvas'
  context : chains
  init : chains.init
  draw : chains.draw
  debug : true

cgl = new PXL.App params
