import gmsh
import sys

'''
This code generates a simple sphere's surface of radius L.
As you can see the way in which gmsh works is a bit cumbersome,
for this reason if you use static environments, it's way easier 
to write the script in python...

NOTE: here no boundary conditions are specificed, the output
is a topologically-closed surface embedded in 3D
'''

# Initialize Gmsh
gmsh.initialize()
gmsh.model.add("sphere")

# Define radius and sphere:
L = 5.0
gmsh.model.occ.addSphere(0,0,0, L)

# Synchronize
gmsh.model.occ.synchronize()

# Set mesh rougnhess (you should test whether normal diffusion is well enough approximated on it)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.05) # Min distance between points of the mesh
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.1) # Max ...

# There are more recent versions, but MFEM currently works only with v2.2
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# Generate mesh
gmsh.model.mesh.generate(2)

# Save mesh to "plane_mesh.msh"
#gmsh.write("../input/mesh/sphere_mesh.msh")
gmsh.write("../input/mesh/sphere_mesh.stl")


# Launch GUI to visualize (optional)
if "-show" in sys.argv:
    gmsh.fltk.run()

# Finalize Gmsh
gmsh.finalize()