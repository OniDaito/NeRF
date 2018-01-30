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
real_ca = []
real_ca.push(new PXL.Math.Vec3(31.89, 53.538, -2.462))
real_ca.push(new PXL.Math.Vec3(29.323, 54.052, -0.956))
real_ca.push(new PXL.Math.Vec3(27.71, 57.258, -2.27))
real_ca.push(new PXL.Math.Vec3(27.985, 57.642, -6.042))

# mod these numbers to match the reference frame

d0 = PXL.Math.Vec3.sub(real_ca[1], real_ca[0])
d1 = PXL.Math.Vec3.sub(real_ca[2], real_ca[1])
d2 = PXL.Math.Vec3.sub(real_ca[3], real_ca[2])

l0 = d0.length()
l1 = d1.length()
l2 = d2.length()

d0.normalize()
d1.normalize()
d2.normalize()

# Starting from zero, here are the real values
final_ca = []
final_ca.push(new PXL.Math.Vec3(0, 0, 0))
final_ca.push(PXL.Math.Vec3.add(final_ca[0], d0.multScalar(l0)))
final_ca.push(PXL.Math.Vec3.add(final_ca[1], d1.multScalar(l1)))
final_ca.push(PXL.Math.Vec3.add(final_ca[2], d2.multScalar(l2)))

# Place holder for computed
computed_ca = []

class TestChain

  _bond_rot : (sp, ep) ->
    dp = PXL.Math.Vec3.sub(ep,sp).normalize()
    y = new PXL.Math.Vec3(0,1,0)
    xp = PXL.Math.Vec3.cross(y,dp)
    dd = Math.acos(dp.dot(y))
    m = new PXL.Math.Matrix4()
    m.rotate(xp, dd)

  constructor : (ca_pos)  ->
    bond_geom = new PXL.Geometry.Cylinder(0.13,50,1,3.82)
    atom_geom = new PXL.Geometry.Sphere(0.5,10)

    pg = new PXL.Colour.RGBA(0.8,0.8,0.8,1.0)
    bond_mat = new PXL.Material.BasicColourMaterial(pg)
    atom_mat = new PXL.Material.BasicColourMaterial(pg)

    @top_node = new PXL.Node()

    residue_atoms_node = new PXL.Node()
    residue_bonds_node = new PXL.Node()

    residue_atoms_node.add(atom_mat)
    residue_bonds_node.add(bond_mat)

    @top_node.add residue_atoms_node
    @top_node.add residue_bonds_node

    idx = 0
    for a in ca_pos
      atom_node = new PXL.Node(atom_geom)
      residue_atoms_node.add(atom_node)
      atom_node.matrix.translate(a)

      if idx != 0
        bond_node = new PXL.Node(bond_geom)
        residue_bonds_node.add(bond_node)
        mp = PXL.Math.Vec3.add(a, ca_pos[idx-1]).multScalar(0.5)
        mm = @_bond_rot(a, ca_pos[idx-1])
        bond_node.matrix.translate(mp).mult(mm)

      idx += 1


class Residue
  constructor : (phi, psi, omega, bond_geom, atom_geom, bond_mat, atom_mat, show_bond) ->
    # Assuming fixed bond lengths and angles with C as the central point
    # Start left handed - initial positions
    @a = new PXL.Math.Vec3(-2.098,1.23,0)
    @b = new PXL.Math.Vec3(-1.33,0,0)
    @c = new PXL.Math.Vec3(0,0,0)

    @phi = phi
    @psi = psi
    @omega = omega

    @c_alpha  = @a
    @nitrogen = @b
    @carbon   = @c

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

    @res_node = new PXL.Node()

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
  next_pos : (prev_res, flip) ->
    a = prev_res.a.clone()
    b = prev_res.b.clone()
    c = prev_res.c.clone()
    d = @a
    na = [@b,@c]
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
        computed_ca.push(d)

    @set_positions()

# Main class for dealing with our 3D chains
class Chains

  _bond_rot : (sp, ep) ->
    dp = PXL.Math.Vec3.sub(ep,sp).normalize()
    y = new PXL.Math.Vec3(0,1,0)
    xp = PXL.Math.Vec3.cross(y,dp)
    dd = Math.acos(dp.dot(y))
    m = new PXL.Math.Matrix4()
    m.rotate(xp, dd)

  _create_chain : (idx) ->
    @residues = []
    bond_geom = new PXL.Geometry.Cylinder(0.13,50,1,3.82)
    atom_geom = new PXL.Geometry.Sphere(0.5,10)

    model_node = new PXL.Node()
    num_residues = @data[idx]['residues'].length

    flip = 1.0
    prev_res = null
    show_bond = false

    for i in [0..num_residues-1]
    #for i in [0..2]
      phi = @data[idx]['angles'][i]['phi']
      psi = @data[idx]['angles'][i]['psi']
      omega = @data[idx]['angles'][i]['omega']

      phi = PXL.Math.degToRad(phi)
      psi = PXL.Math.degToRad(psi)
      omega = PXL.Math.degToRad(omega)

      residue = new Residue(phi, psi, omega, bond_geom, atom_geom, @_get_material_bond(i, num_residues), @_get_material_atom(i, num_residues),show_bond)

      if i > 0
        residue.next_pos(prev_res, flip)

        # Now work on the bonds
        mp = PXL.Math.Vec3.add(residue.a, prev_res.a).multScalar(0.5)
        mm = @_bond_rot(residue.a, prev_res.a)
        residue.bond_node_a.matrix.translate(mp).mult(mm)

      rn = residue.residue_node
      model_node.add rn
      @residues.push residue
      flip *= -1.0

      show_bond = true
      prev_res = residue

    model_node

  _get_material_atom : (i, num_residues) ->
    pink = new PXL.Colour.RGBA(0.8,0.4,0.4,1.0)
    pg = pink.clone()
    pg.r = pg.r / num_residues * (i+1)
    pg.g = pg.g / num_residues * (i+1)
    pg.b = pg.b / num_residues * (i+1)
    calpha_material = new PXL.Material.BasicColourMaterial(pg)

  _get_material_bond : (i, num_residues) ->
    green = new PXL.Colour.RGBA(0.1,0.8,0.1,1.0)
    tg = green.clone()
    tg.r = tg.r / num_residues * (i+1)
    tg.g = tg.g / num_residues * (i+1)
    tg.b = tg.b / num_residues * (i+1)
    backbone_material = new PXL.Material.BasicColourMaterial(tg)

  _parse_cdr : (data) ->
    # Now parse our CDR
    data = eval '(' + data + ')'
    @data = data.data
    @_setup_3d()

  _error : () ->
    # Damn! Error occured
    alert("Error downloading CDR-H3 File")

  _setup_3d : () ->
    # Create the top node and add our camera
    @top_node = new PXL.Node()
    @c = new PXL.Camera.MousePerspCamera new PXL.Math.Vec3(0,0,25)
    @top_node.add @c

    num_models = @data.length
    console.log "num models:" + num_models

    tidx = 0

    # For now just pick the one model
    #while @data[tidx]['name'] != "1NC2_1"
    while @data[tidx]['name'] != "3C6S_2"
      tidx+=1

    for j in [tidx..tidx]
      model_node = @_create_chain(j)
      @top_node.add model_node

    # Add the test chain
    tc = new TestChain(final_ca)
    @top_node.add(tc.top_node)

    uber = new PXL.GL.UberShader @top_node
    @top_node.add uber

    # Print our test values
    console.log("Real CA positions")
    for a in final_ca
      console.log(a)

    console.log("Computed CA Positions")
    for a in computed_ca
      console.log(a)

  init : () ->
    r  = new PXL.Util.Request("data_angles.json")
    r.get ((data) => @_parse_cdr(data)), @_error
    # Basic GL Functions
    GL.enable(GL.CULL_FACE)
    GL.cullFace(GL.BACK)
    GL.enable(GL.DEPTH_TEST)

  draw : () ->
    # Clear and draw our shapes
    GL.clearColor(0.95, 0.95, 0.95, 1.0)
    GL.clear(GL.COLOR_BUFFER_BIT | GL.DEPTH_BUFFER_BIT)
    if @top_node?
      @top_node.draw()

chains = new Chains()

params =
  canvas : 'webgl-canvas'
  context : chains
  init : chains.init
  draw : chains.draw
  debug : true

cgl = new PXL.App params


