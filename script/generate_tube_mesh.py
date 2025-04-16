import gmsh

tube_length = 10
tube_radius = 2

gmsh.initialize()
gmsh.model.add("tube")

tube = gmsh.model.occ.addCylinder(0,0,0,
                                  0,0,tube_length,
                                  tube_radius)
gmsh.model.occ.synchronize()

gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.1)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.2)

gmsh.model.mesh.generate(2)

gmsh.write("../input/mesh/tube_mesh.stl")

gmsh.fltk.run()
gmsh.finalize()


