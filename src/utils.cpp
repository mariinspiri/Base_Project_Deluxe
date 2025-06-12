//
// Created by spiccioni on 2/20/25.
//
#include "../include/utils.h"
#include "../include/Agent.h"
#include "../include/utils.h"

#define NORMAL_SHIFT +0.02
// NOTE: if you want to have the agents "inside" the surface use -0.02

// Euler rotation formula:
geometrycentral::Vector3 rotateVector(geometrycentral::Vector3 point, geometrycentral::Vector3 center, geometrycentral::Vector3 rotation_axis, double rotation_angle) {
    geometrycentral::Vector3 v = point - center;
    geometrycentral::Vector3 k_cross_v = geometrycentral::cross(rotation_axis, v);

    double projection{geometrycentral::dot(v, rotation_axis)};
    geometrycentral::Vector3 rotated_vector;

    for (int d=0;d < 3;d++) {
        rotated_vector[d] = v[d]*cos(rotation_angle) + k_cross_v[d]*sin(rotation_angle)  + rotation_axis[d]*projection*(1-cos(rotation_angle));
    }
    return rotated_vector + center;
}

void utils::saveAgentsPositionToFile(Space* space, std::vector<Agent>& agents, std::string file_name) {
    int num_agents = agents.size();
    int num_of_segments_for_cirlce = 16;
    int total_points = num_agents* (num_of_segments_for_cirlce+1); // center + perimeter
    int total_polygons = num_agents * num_of_segments_for_cirlce;       // each disk subdivided into n_segments
    int total_connectivity_entries = 4 * total_polygons; // each triangle => "3 idx0 idx1 idx2"

    std::ofstream vtk_file(file_name);
    if (!vtk_file) {
        std::cerr << "Cannot open " << file_name << " for writing.\n";
        return;
    }

    // Write VTK header
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "all_agent_disks_with_color\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET POLYDATA\n";


    vtk_file << "POINTS " << total_points << " float\n";


    // First write every vertex's position:
    float d_theta = 2.0f * static_cast<float>(M_PI) / static_cast<float>(num_of_segments_for_cirlce);

    std::vector<geometrycentral::Vector3> all_points;
    all_points.reserve(total_points);

    geometrycentral::Vector3 z_versor{0.0,0.0,1.0};

    for (int i = 0; i < num_agents; i++) {
        const auto& center = agents[i].getGlobalPosition();
        float radius = agents[i].getAgentRadius();

        std::vector<geometrycentral::Vector3> basis = space->gc_getSurfaceTangentBasis(agents[i].gc_position,true);
        geometrycentral::Vector3 normal = basis[2];
        geometrycentral::Vector3 rotation_axis = geometrycentral::cross(z_versor, normal);
        rotation_axis = rotation_axis.normalize();
        double rotation_angle = geometrycentral::angle(z_versor, normal);

        // Center of mass:
        all_points.push_back(center);


        // Perimeter points
        for (int seg = 0; seg < num_of_segments_for_cirlce; seg++) {
            float angle = seg * d_theta;
            float x_off = radius * std::cos(angle);
            float y_off = radius * std::sin(angle);
            geometrycentral::Vector3 perimeter_pos(
                center[0] + x_off,
                center[1] + y_off,
                center[2]
            );

            // now rotate into the surface's tangent plane:
            if (fabs(rotation_angle) > 1e-7 && fabs(rotation_angle - M_PI) > 1e-7) {
                perimeter_pos = rotateVector(perimeter_pos, center, rotation_axis, rotation_angle);
            }

            if (normal.norm() > 1e-7 && normal.isFinite()) {
                perimeter_pos = perimeter_pos + NORMAL_SHIFT*normal;
            }
            all_points.push_back(perimeter_pos);

            if (std::isnan(perimeter_pos[0])) {
                std::cerr<<"ERROR NaN!: agent position"<<center<<" normal:"<<normal<<" rotation angle:"<<rotation_angle<<std::endl;
                return;
            }
        }
    }

    for (const auto& pt : all_points) {
        vtk_file << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
    }

    // polygons
    vtk_file << "POLYGONS " << total_polygons << " " << total_connectivity_entries << "\n";
    for (int i = 0; i < num_agents; i++) {
        int start_index       = i * (num_of_segments_for_cirlce + 1);
        int center_index      = start_index;
        int first_perim_index = start_index + 1;

        for (int seg = 0; seg < num_of_segments_for_cirlce; seg++) {
            int p0 = first_perim_index + seg;
            int p1 = first_perim_index + ((seg + 1) % num_of_segments_for_cirlce);
            vtk_file << "3 " << center_index << " " << p0 << " " << p1 << "\n";
        }
    }


    // Each point has the unique index as a marker
    vtk_file << "CELL_DATA " << total_polygons << "\n";

    // save agent types for visualization
    vtk_file << "SCALARS agent_type int 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < num_agents; i++) {
        int current_type = agents[i].getAgentType();
        for (int seg = 0; seg < num_of_segments_for_cirlce; seg++) {
            vtk_file << current_type << "\n";
        }
    }

    vtk_file << "SCALARS agent_id int 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < num_agents; i++) {
        int current_id = agents[i].getAgentId();
        for (int seg = 0; seg < num_of_segments_for_cirlce; seg++) {
            vtk_file << current_id << "\n";
        }
    }

    vtk_file.close();
}
void utils::saveAgentsFacesToFile(Space* space,std::vector<Agent>& agents,std::string file_name){
    // Open file
    std::ofstream vtk_file(file_name);
    if (!vtk_file) {
        std::cerr << "Cannot open " << file_name << " for writing.\n";
        return;
    }

    auto& mesh     = *(space->gc_mesh);
    auto& geometry = *(space->gc_geometry);

    std::vector<geometrycentral::Vector3> all_positions;
    all_positions.reserve(mesh.nFaces() * 3); // rough guess; each face might have 3 or more vertices

    std::vector<std::vector<int>> polygons;
    polygons.reserve(mesh.nFaces());

    std::vector<int> agent_IDs;    // parallel to polygons
    std::vector<int> agent_types;  // parallel to polygons

    int current_vertex_offset = 0;

    for (auto& ag : agents) {

        int agent_id   = ag.getAgentId();
        int agent_type = ag.getAgentType();

        auto face_set = ag.findFacesWithinRadius();

        for (auto face : face_set) {

            geometrycentral::Vector3 face_normal = geometry.faceNormals[face];
            std::vector<int> face_indices;
            face_indices.reserve(face.degree());

            for (auto he : face.adjacentHalfedges()) {
                auto v = he.vertex();

                geometrycentral::Vector3 pos = geometry.inputVertexPositions[v];

                pos = pos + NORMAL_SHIFT * face_normal;

                all_positions.push_back(pos);

                face_indices.push_back(current_vertex_offset);
                current_vertex_offset++;
            }

            polygons.push_back(face_indices);

            agent_IDs.push_back(agent_id);
            agent_types.push_back(agent_type);
        }
    }

    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "faces_occupied_by_agents\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET POLYDATA\n";

    vtk_file << "POINTS " << all_positions.size() << " float\n";
    for (auto& pt : all_positions) {
        vtk_file << pt[0] << " " << pt[1] << " " << pt[2] << "\n";
    }

    int total_polygons = (int)polygons.size();
    int total_connectivity = 0;
    for (auto& poly : polygons) {
        total_connectivity += (1 + (int)poly.size());
    }

    vtk_file << "POLYGONS " << total_polygons << " " << total_connectivity << "\n";
    for (auto& poly : polygons) {
        vtk_file << poly.size(); // number of vertices in this face
        for (int idx : poly) {
            vtk_file << " " << idx;
        }
        vtk_file << "\n";
    }

    vtk_file << "CELL_DATA " << total_polygons << "\n";

    // Agent IDs
    vtk_file << "SCALARS agent_id int 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (auto& agent_id : agent_IDs) {
        vtk_file << agent_id << "\n";
    }

    vtk_file << "SCALARS agent_type int 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (auto& agent_type : agent_types) {
        vtk_file << agent_type << "\n";
    }
    vtk_file.close();
}

void utils::saveFieldToFile(Space* space, Field& field, std::string file_name, std::string field_name) {
    mfem::GridFunction &f = field.getField();

    std::ofstream output_file(file_name);
    // we can use mfem's built in calls!
    // this one saves the mesh information to file:
    space->mfem_mesh->PrintVTK(output_file, 0, false);

    // this one saves the scalar field information:
    f.SaveVTK(output_file, field_name, 0);

    output_file.close();
}

void utils::saveClearanceTimesToFile(const std::vector<SimulationResult>& results, const std::string &filename)
{
    std::ofstream file(filename);

    if (!file) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    file << "simulation_id, number_best, number_good, clearance_time\n";
    file << std::fixed << std::setprecision(1);
    
    for (size_t i = 0; i < results.size(); ++i) {
        const SimulationResult& res = results[i];
        file << res.simulation_index << "," 
             << res.num_best_cells << ","
             << res.num_good_cells << ","
             << res.clearance_time << "\n";
    }

    file.close();
}
namespace utils {
    void saveTauAndDiffusionCoeffToFile(const std::vector<double>& values_persistence_times, const std::vector<double>& diffCoeff, const std::string& filename) {
        std::ofstream out(filename);
        for (size_t i = 0; i < diffCoeff.size(); ++i) {
            out << values_persistence_times[i] << " " << diffCoeff[i] << "\n";
        }
        out.close();
    }
}