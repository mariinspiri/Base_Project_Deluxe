#include <algorithm>
#include <chrono>
#include <iostream>
#include <geometrycentral/surface/manifold_surface_mesh.h>

#include "include/utils.h"

#include "include/Agent.h"
#include "include/CollisionManager.h"
#include "include/Field.h"
#include "include/Space.h"

// NOTE-CONVENTION: good cells are positive, bad one ar negative as obvious
#define AGENT_TYPE_GOOD_CELL +1
#define AGENT_TYPE_EVIL_CELL -1

using namespace std;
using namespace geometrycentral::surface;

int main(int argc, char *argv[]) {
    auto time_start = std::chrono::high_resolution_clock::now();

    // load mesh and declare space
    string mesh_file("plane_mesh");
    Space space("../input/mesh/"+mesh_file+".stl", 91);

    int n_of_elements = space.gc_mesh->nFaces();

    // uniform distribution [0, 1]
    std::uniform_real_distribution<double> unif_dist(0.0, 1);

    vector<Agent> agents;
    // Agents initialization:
    {
        // agents will contain every type of agent
        int agent_id_counter {1};

        // initialize good cells:
        int num_good_cells {10};
        double persistence_time {2};
        double radius_good_cells {0.2};
        double speed_good_cells {0.2};

        for (int i = 0;i < num_good_cells;i++) {
            // random face, in the center of the triangle:
            Face face = space.gc_mesh->face(static_cast<int>(n_of_elements * unif_dist(space.gen)));
            SurfacePoint initial_position(face, {1./3., 1./3., 1./3.});

            Agent agent(&space, initial_position, radius_good_cells,
                                            agent_id_counter,
                                            AGENT_TYPE_GOOD_CELL,
                                            persistence_time,
                                            persistence_time*unif_dist(space.gen)); //initial phase-lag
            agent.setLocalVelocity({speed_good_cells, 0});
            agents.push_back(agent);
            ++agent_id_counter;
        }

        // initialize "evil" cells
        int num_evil_cell  {10};
        double evil_radius {0.05};

        for (int i = 0; i < num_evil_cell; i++) {
            Face face = space.gc_mesh->face(static_cast<int>(n_of_elements*unif_dist(space.gen)));
            SurfacePoint initial_position(face, {1./3., 1./3., 1./3.});

            Agent agent(&space, initial_position, evil_radius,
                                            agent_id_counter,
                                            AGENT_TYPE_EVIL_CELL);
            agents.push_back(agent);
            ++agent_id_counter;
        }
    }

    // checks that there are no overlaps... (NOTE: if some evil one are already touching the good ones, they could be killed)
    CollisionManager collision_manager(&space);
    {
        collision_manager.checkCollisions(agents);
        collision_manager.fixCollisions(agents);
    }

    double T {60};          // minutes
    double dt {0.1};        // minutes

    double D {0.01};        // Diffusion coefficient
    double lambda  {0.01};   // Reaction/decay rate
    double S_0 {1};         // source intensity (NOTE: in the PDE it should scale with the area e.g. (S_0/Area), so it would be like having S_0=|E|^-1)

    // initialize the chemokines' scalar field
    Field chemokines(&space);
    chemokines.setup(dt, D, lambda);
    chemokines.computeSources(agents, S_0);

    int file_index {0};

    for (int i = 0; i < static_cast<int>(T/dt); i++) {
        cout<<" step:"<<i<<" time:"<<(i*dt)<<endl;

        // shuffle the agents:
        std::shuffle(agents.begin(), agents.end(), space.gen);

        // AGENTS DYNAMIC:
        for (auto &agent : agents) {
            // NOTE: I chose to put the calls directly here, but you can think of implementing a agent.doStep(dt) function to make the main cleaner
            if (agent.getAgentType() == AGENT_TYPE_GOOD_CELL) {
                agent.move(dt);

                // calling the persistence timer updates the internal one - returns true and reset itself when it reaches 0
                if (agent.persistenceTimer(dt)){
                    agent.computeNewBPRWVelocity();
                }
            }
        }

        // CHECK "COLLISIONS" (INTERACTIONS OF ALL SORTS):
        collision_manager.checkCollisions(agents);
        bool agents_list_updated = collision_manager.fixCollisions(agents);

        // UPDATE THE CHEMOKINES:
        // update the source term if some evil cells were actually removed
        if (agents_list_updated) {
            chemokines.computeSources(agents, S_0);
        }

        // update the sinks and save inside the amount of bound receptors inside of the good agents (along with the current average gradient)
        chemokines.computeSinksAndBindReceptors(dt, agents);
        chemokines.step(dt);

        // visualizations:
        if (i % 1 == 0){
            utils::saveAgentsPositionToFile(&space, agents,
                                            "../output/" + mesh_file + "_agents_disks_step_" + std::to_string(file_index) + ".vtk");

            //utils::saveAgentsFacesToFile(&space, agents, "../output/" + mesh_file + "_agents_faces_step_" + std::to_string(file_index) + ".vtk");

            utils::saveFieldToFile(&space, chemokines,
                                    "../output/"+ mesh_file + "_field_step_"+ std::to_string(file_index)+ ".vtk", "chemokines");

            ++file_index;
        }
    }

    auto time_end = std::chrono::high_resolution_clock::now();
    auto time_duration = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count()/1000;

    // Print the execution time
    std::cout << "Execution time: " << time_duration << " s" << std::endl;
    return 0;
}