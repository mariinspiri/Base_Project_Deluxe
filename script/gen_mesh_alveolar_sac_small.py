import gmsh
import math
import numpy as np
from numpy.linalg import norm

# Parameters
tube_radius = 0.5        # Radius of the airway tube
tube_length = 2        # Length of the tube
alveolus_radius = 1.5
overlap_factor = 0.8
overlap_factor_lateral=1.1
pore_radius = 0.3
wall_thickness = 0.05


def aggiungi_acino(new_pos, positions, cumulative_result, wall_thickness=wall_thickness, acino_radius=alveolus_radius,pore_radius=0):
    new_acino = gmsh.model.occ.addSphere(new_pos[0],
                                         new_pos[1],
                                         new_pos[2], acino_radius)
    acino_tag = new_acino
    selected_regions_collector = []

    for pos in positions:
        if norm(pos - new_pos) < 2*acino_radius and norm(pos - new_pos) > 0:
            mid_pos = (new_pos + pos)/2
            directional_vector = new_pos - pos
            directional_vector = directional_vector/norm(directional_vector)

            size = 1500
            box = gmsh.model.occ.addBox(-size,-size,0,
                                        2*size,2*size,2*size)
            
            complementary_box = gmsh.model.occ.addBox(-size,-size,0,
                                        2*size,2*size,-2*size)

            rotational_axis = np.array([-directional_vector[1], directional_vector[0],0])
            angle = np.arccos(np.dot(np.array([0,0,1]), directional_vector))

            if norm(rotational_axis) > 10e-6:
                gmsh.model.occ.rotate([(3,box)], 0,0,0, 
                                    rotational_axis[0], rotational_axis[1], rotational_axis[2],
                                    angle)
                gmsh.model.occ.rotate([(3,complementary_box)], 0,0,0, 
                                    rotational_axis[0], rotational_axis[1], rotational_axis[2],
                                    angle)
            
            translation_vector =  mid_pos + directional_vector*(wall_thickness/2)
            gmsh.model.occ.translate([(3,box)], translation_vector[0], translation_vector[1], translation_vector[2])

            translation_vector_compl =  mid_pos - directional_vector*(wall_thickness/2)
            gmsh.model.occ.translate([(3,complementary_box)], translation_vector_compl[0], translation_vector_compl[1], translation_vector_compl[2])

            new_acino,_ = gmsh.model.occ.intersect([(3,acino_tag)], [(3,box)])
            acino_tag = new_acino[0][1]
            
            selected_regions_collector.append((3, complementary_box))
            gmsh.model.occ.synchronize()

    good_region = [selected_regions_collector[0]]
    if len(selected_regions_collector) > 1:
        good_region,_ = gmsh.model.occ.fuse([selected_regions_collector[0]], selected_regions_collector[1:])

    cumulative_result,_ = gmsh.model.occ.intersect(cumulative_result, good_region)


    #now put the holes:
    for pos in positions:
        if norm(pos - new_pos) < 2*acino_radius and norm(pos - new_pos) > 0:
            mid_pos = (new_pos + pos)/2
            directional_vector = new_pos - pos
            directional_vector = directional_vector/norm(directional_vector)
            rotational_axis = np.array([-directional_vector[1], directional_vector[0],0])
            angle = np.arccos(np.dot(np.array([0,0,1]), directional_vector))

            cylinder_start = mid_pos - (wall_thickness/2)*directional_vector
            connecting_cylinder  = gmsh.model.occ.addCylinder(cylinder_start[0],
                                                            cylinder_start[1],
                                                            cylinder_start[2],
                                                            0,0,
                                                            wall_thickness,pore_radius
                                                            )
            
            rotational_axis = np.array([-directional_vector[1], directional_vector[0],0])
            angle = np.arccos(np.dot(np.array([0,0,1]), directional_vector))

            if norm(rotational_axis) > 10e-6:
                gmsh.model.occ.rotate([(3,connecting_cylinder)], cylinder_start[0], cylinder_start[1],cylinder_start[2], 
                                        rotational_axis[0], rotational_axis[1], rotational_axis[2],
                                        angle)
            new_acino,_ = gmsh.model.occ.fragment(new_acino, [(3,connecting_cylinder)])
        
    cumulative_result,_ = gmsh.model.occ.fuse(cumulative_result, new_acino)
    return cumulative_result


def generate_bottomup_model():
    gmsh.initialize()
    gmsh.model.add("alveolar_bottomup_v3")

    #tube
    tube = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, tube_length,tube_radius)
    gmsh.model.occ.synchronize()


    #position defs.
    z_shift = math.sqrt(alveolus_radius**2 - tube_radius**2) #such that it's like putting a sphere on a cylinder

    pos_central = np.array([0, 0, tube_length+z_shift])
    pos_bottom = np.array([0, 0, tube_length+z_shift+2*overlap_factor*alveolus_radius])
    pos_right = np.array([overlap_factor_lateral*alveolus_radius, 0, tube_length+z_shift+overlap_factor*alveolus_radius])
    pos_left = np.array([-overlap_factor_lateral*alveolus_radius, 0, tube_length+z_shift+overlap_factor*alveolus_radius])
    # pos_in = np.array([0, overlap_factor_lateral*alveolus_radius, tube_length+z_shift+overlap_factor*alveolus_radius])
    # pos_out = np.array([0, -overlap_factor_lateral*alveolus_radius, tube_length+z_shift+overlap_factor*alveolus_radius])

    sphere_centers = [pos_central, pos_bottom, pos_right, pos_left] #, pos_in, pos_out]

    central_sphere = gmsh.model.occ.addSphere(pos_central[0], 
                                              pos_central[1],
                                              pos_central[2], alveolus_radius)
    gmsh.model.occ.synchronize()


    cumulative_result,_ = gmsh.model.occ.fuse([(3,tube)], [(3,central_sphere)])
    gmsh.model.occ.synchronize()
    
    for k in range(1,len(sphere_centers)):
        cumulative_result = aggiungi_acino(sphere_centers[k], sphere_centers, cumulative_result, pore_radius=pore_radius)
        gmsh.model.occ.synchronize()
    
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.05)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

    gmsh.model.mesh.generate(2)


    #gmsh.write("../input/mesh/alveolar_sac_mesh_small.msh")
    gmsh.write("../input/mesh/alveolar_sac_mesh_small.stl")
    gmsh.finalize()

if __name__ == "__main__":
    generate_bottomup_model()
