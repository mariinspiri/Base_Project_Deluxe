//
// Created by spiccioni on 2/20/25.
//
#include "../include/Space.h"

using namespace mfem;
using namespace std;

using namespace geometrycentral;
using namespace geometrycentral::surface;

Space::Space(std::string gc_mesh_file_name) {
    loadGCMeshFromFile(gc_mesh_file_name);
    createMFEMMeshFromGC();
}

Space::~Space()
{
    delete finite_element_space;
    delete finite_element_collection;
    delete mfem_mesh;
}


void Space::loadGCMeshFromFile(std::string filename) {
    std::tie(gc_mesh, gc_geometry) = readManifoldSurfaceMesh(filename);

    // now initialize geometrical information about the mesh:
    gc_geometry->requireFaceNormals();
    gc_geometry->requireEdgeLengths();
    gc_geometry->requireVertexNormals();
}

void Space::createMFEMMeshFromGC(){
    // GC-mesh -> MFEM:
    // this will also ensure that they share the same indexes for the vertices

    int n_vertexes = gc_mesh->nVertices();
    int n_faces = gc_mesh->nFaces();

    static std::vector<double> coords;
    coords.resize(3 * n_vertexes);

    for (geometrycentral::surface::Vertex v : gc_mesh->vertices()) {
        size_t i = v.getIndex();
        auto p = gc_geometry->inputVertexPositions[v];
        coords[3*i+0] = p.x;
        coords[3*i+1] = p.y;
        coords[3*i+2] = p.z;
    }

    std::vector<int> face_vertices(3 * n_faces);
    for (Face f : gc_mesh->faces()) {
        size_t face_idx = f.getIndex();
        int k = 0;
        for (geometrycentral::surface::Vertex fv : f.adjacentVertices()) {
            face_vertices[3*face_idx + k] = fv.getIndex();
            k++;
        }
    }

    // declare the mesh
    mfem_mesh = new Mesh(manifold_dimension, n_vertexes, n_faces, 0, space_dimension);

    // add nodes:
    for (int i = 0; i < n_vertexes; i++) {
        mfem_mesh->AddVertex(coords[3*i + 0],
                           coords[3*i + 1],
                           coords[3*i + 2]);
    }
    // add the faces:
    for (int f = 0; f < n_faces; f++) {
        int v0 = face_vertices[3*f + 0];
        int v1 = face_vertices[3*f + 1];
        int v2 = face_vertices[3*f + 2];
        mfem_mesh->AddTriangle(v0, v1, v2);
    }

    mfem_mesh->FinalizeTopology();

    // initialize the functional space:
    finite_element_collection = new mfem::H1_FECollection(functions_order, manifold_dimension);
    finite_element_space = new mfem::FiniteElementSpace(mfem_mesh, finite_element_collection);
}

std::vector<Vector3> Space::gc_getSurfaceTangentBasis(SurfacePoint surface_point, bool include_normal=false) {
    Vector3 e1, e2, normal; // Local tangent space basis

    // Vertex case:
    // if you are on a vertex we can use the adjacent edges to get the tangent plane
    // (it's better than just picking one of the faces)
    if (surface_point.type == SurfacePointType::Vertex) {
        geometrycentral::surface::Vertex v = surface_point.vertex;
        normal = gc_geometry->vertexNormals[v];

        // Compute a local tangent basis using adjacent edges
        for (Halfedge he : v.outgoingHalfedges()) {
            Vector3 candidate = (gc_geometry->inputVertexPositions[he.twin().vertex()] -
                                 gc_geometry->inputVertexPositions[v]).normalize();
            if (fabs(dot(candidate, normal)) < 1e-6) {  // Ensure perpendicularity
                e1 = candidate;
                break;
            }
        }
        e2 = cross(normal, e1).normalize();

    }
    // Edge case:
    // if you are on an edge we use one of the faces connected to it to get the normal
    else if (surface_point.type == SurfacePointType::Edge) {
        Edge e = surface_point.edge;
        Halfedge he = e.halfedge();
        Vector3 v0 = gc_geometry->inputVertexPositions[he.vertex()];
        Vector3 v1 = gc_geometry->inputVertexPositions[he.twin().vertex()];

        e1 = (v1 - v0).normalize();  // Tangent along the edge

        // Compute normal from one of the adjacent faces
        Face f = he.face();
        normal = f.isBoundaryLoop() ? Vector3{0, 0, 0} : gc_geometry->faceNormal(f);

        e2 = cross(normal, e1).normalize();  // Perpendicular to the edge

    }
    // Face case:
    // here we construct the basis from the face's vertices
    else if (surface_point.type == SurfacePointType::Face) {
        Face f = surface_point.face;
        Halfedge he = f.halfedge();
        Vector3 v0 = gc_geometry->inputVertexPositions[he.vertex()];
        Vector3 v1 = gc_geometry->inputVertexPositions[he.next().vertex()];

        e1 = (v1 - v0).normalize();
        normal = gc_geometry->faceNormal(f);
        e2 = cross(normal, e1).normalize();
    } else {
        throw std::runtime_error("Invalid SurfacePoint type.");
    }

    if (include_normal) {
        return {e1,e2,normal};
    }

    return {e1, e2};
}


void Space::convertGlobalVectorToLocalVector(SurfacePoint s_point, Vector3 global_vec, Vector2 &local_v) {
    std::vector<Vector3> tangent_basis = gc_getSurfaceTangentBasis(s_point);

    // Project the global vector onto the local basis
    double localX = dot(global_vec, tangent_basis[0]);
    double localY = dot(global_vec, tangent_basis[1]);

    local_v.x = localX;
    local_v.y = localY;
}

void Space::convertLocalVectorToGlobalVector(SurfacePoint s_point, Vector2 local_vec, Vector3 &global_vec) {
    std::vector<Vector3> tangent_basis = gc_getSurfaceTangentBasis(s_point);

    global_vec = local_vec.x*tangent_basis[0] + local_vec.y*tangent_basis[1];
}
