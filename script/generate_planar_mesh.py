import gmsh
import sys

'''
This code generates a simple mesh of a plane (L x L).
As you can see the way in which gmsh works is a bit cumbersome,
for this reason if you use static environments, it's way easier 
to write the script in python...

NOTE: here no boundary conditions are specificed, the output
is a topologically-closed surface embedded in 2D
'''

# Initialize Gmsh
gmsh.initialize()
gmsh.model.add("plane_mesh")

# Define corner points of the unit square
L = 10.0  # Length of square
p1 = gmsh.model.geo.addPoint(0, 0, 0)
p2 = gmsh.model.geo.addPoint(L, 0, 0)
p3 = gmsh.model.geo.addPoint(L, L, 0)
p4 = gmsh.model.geo.addPoint(0, L, 0)

# Define edges
l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)

# Define curve loop and surface
cl = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
s = gmsh.model.geo.addPlaneSurface([cl])

# Synchronize
gmsh.model.geo.synchronize()

# Set mesh rougnhess (you should test whether normal diffusion is well enough approximated on it)
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.05) # Min distance between points of the mesh
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.1) # Max ...

# There are more recent versions, but MFEM currently works only with v2.2
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# Generate mesh
gmsh.model.mesh.generate(2)

# Save mesh to "plane_mesh.msh"
#gmsh.write("../input/mesh/plane_mesh.msh")
gmsh.write("../input/mesh/plane_mesh.stl")

# Launch GUI to visualize (optional)
if "-show" in sys.argv:
    gmsh.fltk.run()

# Finalize Gmsh
gmsh.finalize()