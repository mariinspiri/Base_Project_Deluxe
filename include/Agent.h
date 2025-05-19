//
// Created by spiccioni on 2/20/25.
//

#ifndef AGENT_H
#define AGENT_H

#include <geometrycentral/surface/trace_geodesic.h>

#include "Space.h"


class Agent {
public:
    Agent(Space *s, geometrycentral::surface::SurfacePoint &initial_position,
            double radius, int agent_id, int agent_type=0, double persistence_period=1e15,
            double initial_persistence_timer=0.0);

    // GET methods:
    int getAgentId(){return agent_id;}
    int getAgentType(){return agent_type;}
    double getAgentRadius(){return radius;}
    geometrycentral::Vector2 getLocalVelocity();
    geometrycentral::Vector3 getGlobalVelocity();
    geometrycentral::Vector3 getGlobalPosition();

    // move the agent along geodesics, if use_speed uses also the previously saved speed, otherwise sets to 1
    int move(double dt, bool use_speed=true, int recursion_index=0);

    // updates the internal time, returns true if timer reached 0
    bool persistenceTimer(double dt);

    // randomly rotate the internal velocity or follow the gradient
    void computeNewBPRWVelocity();

    // converts a global velocity into local one and stores the value
    void setLocalVelocity(std::vector<double> velocity);

    // GC variables for position and velocity:
    geometrycentral::surface::SurfacePoint gc_position;

    geometrycentral::Vector2 gc_local_velocity_direction{0., 0.};
    double speed = 0.0;

    // find faces within a certain radius, if nothing is specified the internal radius is used:
    std::unordered_set<geometrycentral::surface::Face> findFacesWithinRadius(double coverage_radius=-1);

    // Ligand receptors dynamic:
    double getFreeReceptors(){return R_concentration;}
    void stepLigandReceptors(double dt);
    void updateLigandReceptors(double bound_receptors, geometrycentral::Vector3 average_gradient);

    // parameters:
    double k_binding{0.01};
    double k_internalized{0.07};
    double k_recycled{0.05};
    double sensitivity_to_gradient{5*1e2}; // what in Pollmächer is sigma

private:
    // Pointer to the belonging space:
    Space* space;

    // Agent information:
    int agent_type{0};
    int agent_id{-1};

    // Physical properties:
    double radius{0.0};

    // PRW period and internal timer
    double persistence_period{1e15};
    double persistence_timer{1.0};

    // Ligand-Receptor Dynamics:
    double total_concentration_sum{1e4};

    double R_concentration{total_concentration_sum};
    double LR_concentration{0};
    double R_star_concentration{0};

    double newly_bound_receptors{0.0};
    geometrycentral::Vector3 newly_average_gradient{0., 0., 0.};

    geometrycentral::Vector3 cumulative_gradient{0., 0., 0.};

    // GC option helper:
    geometrycentral::surface::TraceOptions options = geometrycentral::surface::defaultTraceOptions;

    // memory for covered faces:
    std::unordered_set<geometrycentral::surface::Face> memory_covered_faces;
    bool updated_covered_faces = false;
};
#endif //AGENT_H
