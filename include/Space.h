//
// Created by spiccioni on 2/20/25.
//
#ifndef SPACE_H
#define SPACE_H

#include <iostream>
#include <mfem.hpp>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <geometrycentral/surface/vertex_position_geometry.h>
#include <geometrycentral/surface/meshio.h>
#include <geometrycentral/surface/surface_point.h>

// This class is essentially a container for all the shared data structures: GC's mesh and geometry, MFEM's mesh and functional space...
class Space {
public:
    // this constructor uploads the mesh file into GC's mesh and MFEM's one + declares the unique random generator used through out all of the simulation
    Space(std::string gc_mesh_file_name, int seed);
    ~Space();

    // loads mesh from file ".stl" format, and declares gc_mesh/gc_geometry
    void loadGCMeshFromFile(std::string filename);

    // from gc_mesh and gc_geometry constructs MFEM's mesh
    void createMFEMMeshFromGC();

    // FUNCTIONS to work with local vs global tangent spaces:
    // calculates a basis for T_{surface_point}[MESH]
    std::vector<geometrycentral::Vector3> gc_getSurfaceTangentBasis(geometrycentral::surface::SurfacePoint surface_point, bool include_normal);
    //convert a vector from the T_{surface_point}[MESH] basis to the global one
    void convertLocalVectorToGlobalVector(geometrycentral::surface::SurfacePoint s_point, geometrycentral::Vector2 local_vec, geometrycentral::Vector3 &global_vec);
    // viceversa
    void convertGlobalVectorToLocalVector(geometrycentral::surface::SurfacePoint sPoint,geometrycentral::Vector3 globalVec, geometrycentral::Vector2 &local_v);

    //MESHES:
    // geometry-central (GC) mesh:
    std::unique_ptr<geometrycentral::surface::ManifoldSurfaceMesh> gc_mesh;
    std::unique_ptr<geometrycentral::surface::VertexPositionGeometry> gc_geometry;

    // MFEM mesh:
    mfem::Mesh* mfem_mesh;

    // functional space used to do FEM on the fields
    mfem::FiniteElementCollection* finite_element_collection;    // Type of bump functions (H1, DG, etc.)
    mfem::FiniteElementSpace* finite_element_space;

    std::mt19937 gen;

private:
    // dimensions
    int manifold_dimension {2};
    int space_dimension {3};

    // order of the test functions (1: linear, 2: quadratic)
    int functions_order {1};
};

#endif //SPACE_H
