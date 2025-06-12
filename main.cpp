#include <algorithm>
#include <chrono>
#include <iostream>
#include <geometrycentral/surface/manifold_surface_mesh.h>
#include <utility>
#include <cstdlib> 

#include "include/utils.h"

#include "include/Agent.h"
#include "include/CollisionManager.h"
#include "include/Field.h"
#include "include/Space.h"
#include "include/AgentTypes.h"

using namespace std;
using namespace geometrycentral::surface;

//compute clearance time from simulation
float calculateClearanceTime(int num_best_cells, int num_good_cells, int id_sim){
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

        // initialize the best cells
        
        double persistence_best_time {2};
        double radius_best_cells {0.15};
        double speed_best_cells {0.8};

        for (int i = 0;i < num_best_cells;i++) {
            // random face, in the center of the triangle:
            Face face = space.gc_mesh->face(static_cast<int>(n_of_elements * unif_dist(space.gen)));
            SurfacePoint initial_position(face, {1./3., 1./3., 1./3.});

            Agent agent(&space, initial_position, radius_best_cells,
                                            agent_id_counter,
                                            AGENT_TYPE_BEST_CELL,
                                            persistence_best_time,
                                            persistence_best_time*unif_dist(space.gen),
                                            true); //initial phase-lag
            agent.setLocalVelocity({speed_best_cells, 0});
            agents.push_back(agent);
            ++agent_id_counter;
        }
        
        // initialize good cells:
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
                                            persistence_time*unif_dist(space.gen),
                                            true); //initial phase-lag
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

    double T {20};          // minutes
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
        // shuffle the agents:
        std::shuffle(agents.begin(), agents.end(), space.gen);

        // AGENTS DYNAMIC:
        for (auto &agent : agents) {
            agent.doStep(dt);
        }

        // CHECK "COLLISIONS" (INTERACTIONS OF ALL SORTS):
        collision_manager.checkCollisions(agents);
        bool agents_list_updated = collision_manager.fixCollisions(agents);

        //chack if clearance time was reached
        int evil_cell_count = 0;
        for(int j = 0; j < agents.size(); j++){
            if(agents[j].getAgentType() == AGENT_TYPE_EVIL_CELL){
                evil_cell_count++;
            }
        }
        if(evil_cell_count <= 0){
            return i*dt;
        }

        cout<<" time:"<< (i*dt) <<" evil cells: "<< evil_cell_count <<endl;

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

    //if clearance time was not reached
        return -1;
}

//compute diffusion coefficient from MSD values
double computeDiffusionCoefFromMSD(const std::vector<float>& msd, float dt) {
    int n = msd.size();
    if (n == 0) {
        throw std::invalid_argument("MSD vector is empty");
    }

    double sum_t = 0.0;
    double sum_msd = 0.0;
    double sum_tt = 0.0;
    double sum_tmsd = 0.0;

    for (int i = 0; i < n; ++i) {
        double t = i * dt;
        sum_t += t;
        sum_msd += msd[i];
        sum_tt += t * t;
        sum_tmsd += t * msd[i];
    }

    double denominator = n * sum_tt - sum_t * sum_t;
    if (denominator == 0.0) {
        throw std::runtime_error("Denominator in linear regression calculation is zero");
    }

    double slope = (n * sum_tmsd - sum_t * sum_msd) / denominator;
    double intercept = (sum_msd - slope * sum_t) / n;

    double D = slope / 4.0;//difusion coefficient

    return D;
}

//compute MSD from simulation
double calculateDiffusionCoef(int num_good_cells, double persistence_time, double speed_good_cells){
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
        //double persistence_time {2};
        double radius_good_cells {0.2};
        //double speed_good_cells {0.2};

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

    int file_index {0};

    //Mean Squared Displacement
    std::vector<float> msd;
    //initial position for all agents
    std::vector<std::pair<float, float>> initial_pos;
    initial_pos.resize(num_good_cells);
    for (auto &agent : agents) {
        int id = agent.getAgentId() - 1;
        geometrycentral::Vector3 pos = agent.getGlobalPosition();
        initial_pos[id] = {pos.x, pos.y};
    }

    for (int i = 0; i < static_cast<int>(T/dt); i++) {
        // shuffle the agents:
        std::shuffle(agents.begin(), agents.end(), space.gen);
        
        // AGENTS DYNAMIC:
        for (auto &agent : agents) {
            agent.doStep(dt);
        }

        // CHECK "COLLISIONS" (INTERACTIONS OF ALL SORTS):
        collision_manager.checkCollisions(agents);
        bool agents_list_updated = collision_manager.fixCollisions(agents);

        //calculate msd based on global position
        std::vector<float> local_msd;
        local_msd.resize(num_good_cells);
        for (auto &agent : agents) {
            int id = agent.getAgentId() - 1;
            geometrycentral::Vector3 position = agent.getGlobalPosition();
            
            float dx = position.x - initial_pos[id].first;
            float dy = position.y - initial_pos[id].second;

            local_msd[id]= dx*dx + dy*dy;
        }
        

        float sum = 0.0f;
        for (float val : local_msd) {
            sum += val;
        }
        if (!local_msd.empty()) {
            msd.push_back(sum / local_msd.size());
        }

        // visualizations:
        if (i % 1 == 0){
            utils::saveAgentsPositionToFile(&space, agents,
                                            "../output/" + mesh_file + "_agents_disks_step_" + std::to_string(file_index) + ".vtk");
            
            ++file_index;
        }
        
    }

    auto time_end = std::chrono::high_resolution_clock::now();
    auto time_duration = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_start).count()/1000;

    // Print the execution time
    //std::cout << "Execution time: " << time_duration << " s" << std::endl;

    double diffCoeff = computeDiffusionCoefFromMSD(msd, dt);
    return diffCoeff;
}

int main(int argc, char *argv[]) {
    /*    
    std::vector<double> values_persistence_time = {0.1, 0.5, 1.0, 2.0, 5.0, 10.0};
    std::vector<double> values_speed = {0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0};
    double fixed_speed = 0.2;
    double fixed_tau = 2.0;

    // Vary persistence time, fixed speed
    std::ofstream out1("../output/diffusion_vs_tau.txt");
    out1 << "tau D\n";
    for (double tau : values_persistence_time) {
        double D = calculateDiffusionCoef(10, tau, fixed_speed);
        out1 << tau << " " << D << "\n";
    }
    out1.close();

    // Vary speed, fixed persistence time
    std::ofstream out2("../output/diffusion_vs_speed.txt");
    out2 << "speed D\n";
    for (double v : values_speed) {
        double D = calculateDiffusionCoef(10, fixed_tau, v);
        out2 << v << " " << D << "\n";
    }
    out2.close();

    // Run the plot script
    int plot = system("python3 ../script/plot_diffCoeff.py");
    if (plot != 0) {
        std::cerr << "Error: Python script failed to run.\n";
    }
    */

    
    //num_cells consist of num_best and num_good
    std::vector<std::pair<int,int>> num_cells = {
        //{0, 10},
        //{2, 8},
        {5, 5},
        //{8, 2},
        //{10, 0}
    };

    std::vector<SimulationResult> results;
    int id_sim {0};

    //run simulation several times and get clearance time as a return parameter
    for (const auto& num : num_cells) {
        for(int i = 0; i < 1; i++){
            float clearance_time = calculateClearanceTime(num.first, num.second, id_sim);
            results.push_back({id_sim, num.first, num.second, clearance_time});
            id_sim++;
        }
    }
    //save clearance time to the .csv file and save plotted diagram in output 
    utils::saveClearanceTimesToFile(results, "../output/clearance_times.csv");
    int plot = system("python3 ../script/plot_stats.py");
    if (plot != 0){
        std::cerr<< "Error: Python script failed to run.\n";
    }
    
    return 0;
}