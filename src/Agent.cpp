//
// Created by spiccioni on 2/20/25.
//
#include <queue>

#include "../include/Agent.h"

using namespace std;
using namespace geometrycentral;
using namespace geometrycentral::surface;

const int MAX_NUMBER_OF_MESH_SEARCH_ITERS = 1e5;
struct VertexHash {
    std::size_t operator()(const geometrycentral::surface::Vertex& v) const {
        return std::hash<size_t>()(v.getIndex());
    }
};


Agent::Agent(Space *s, SurfacePoint &initial_position, double agent_radius, int agent_id, int agent_type, double persistence_period, double initial_persistence_timer):
    space(s) {

    this->persistence_period = persistence_period;
    persistence_timer = initial_persistence_timer;
    radius = agent_radius;

    // safeguard for the index:
    if (agent_id == 0) {
        cerr<<"Cannot initialize agent's index to 0!"<<endl;
        return;
    }
    this->agent_id = agent_id;
    this->agent_type = agent_type;

    gc_position = initial_position;
}

int Agent::move(double dt, bool use_speed, int recursion_index) {
    // safeguard against infinite-loops:
    if (recursion_index > 1e2) {
        return -1;
        // or implement resolution strategy
    }

    // if we move we have to re-run the faces-finder function:
    updated_covered_faces = false;

    double gc_speed = speed;
    if (!use_speed){gc_speed=1;}

    TraceGeodesicResult result =traceGeodesic(*space->gc_geometry, gc_position, gc_local_velocity_direction*gc_speed*dt, options);
    gc_position = result.endPoint;
    gc_local_velocity_direction = result.endingDir.normalize();

    double remaining_length = min(gc_speed*dt - result.length, gc_speed*dt);

    //TODO: if you have boundary conditions you can put them here:
    if (result.hitBoundary || remaining_length > 1e-6) {
        Vector3 global_position = getGlobalPosition();


        // save the velocity direction before the corrective-shift:
        Vector3 global_vel_direction;
        space->convertLocalVectorToGlobalVector(gc_position, gc_local_velocity_direction, global_vel_direction);


        {
            // BC: x-direction periodic boundary conditions:
            // NOTE: this is a quick and dirty solution that doesn't require much coding
            //       if your mesh has some more complicated structures, you can think of creating a table containing element-to-element
            //       pairings to encode the periodicity, otherwise you could also perform the shift in global coordinates and then loop
            //       over the elements to check on which one we are, although this might be more efficient if done with MFEM...
            /*
            double periodic_shift=0.0;
            if (global_position[0] >= 10.0) {
                periodic_shift = -10.0;
            }else if (global_position[0] <= 0.0) {
                periodic_shift = 10.0;
            }
            if (periodic_shift != 0.0) {
                // you have to find the new element and the new local coordinates:
                Vector3 shift{periodic_shift,0,0};
                Vector2 local_shift;

                space->convertGlobalToLocalVector(gc_position, shift, local_shift);
                result = traceGeodesic(*space->gc_geometry, gc_position, local_shift, options);
                gc_position = result.endPoint;

                space->convertGlobalToLocalVector(gc_position, global_vel_direction, gc_local_velocity_direction);
            }
            */
        }

        // REFLECTIVE BOUNDARY CONDITIONS:
        // NOTE: I use an offset of 1e-12 to avoid the interpolation errors introduced by getGlobalPosition()
        // BC: x-direction reflective boundary conditions:
        if (10.0 - global_position[0] <= 1e-12 || global_position[0] <= 1e-12) {
            global_vel_direction[0] *= -1.0;
        }

        // BC: y-direction reflective boundary conditions:
        if (10.0  - global_position[1] <= 1e-12 || global_position[1] <= 1e-12) {
            global_vel_direction[1] *= -1.0;
        }

        // now update the velocity:
        space->convertGlobalVectorToLocalVector(gc_position, global_vel_direction, gc_local_velocity_direction);
        gc_local_velocity_direction = gc_local_velocity_direction.normalize();

        // move for how much is left:
        return move(remaining_length, false, recursion_index+1);
    }
    return 0;
}


bool Agent::persistenceTimer(double dt) {
    persistence_timer -= dt;
    if (persistence_timer <= 0) {
        persistence_timer = persistence_period;
        return true;
    }
    return false;
}

std::unordered_set<Face> Agent::findFacesWithinRadius(double coverage_radius) {
    if (coverage_radius == -1) {
        coverage_radius = radius;

        if (updated_covered_faces) {
            return memory_covered_faces;
        }
    }
    std::unordered_set<Face> faces_within_radius;

    // Priority queue for Dijkstra's algorithm
    std::priority_queue<std::pair<double, Vertex>, std::vector<std::pair<double, Vertex>>, std::greater<>> priority_queue;

    // Distance map with custom hash function
    std::unordered_map<Vertex, double, VertexHash> distance;

    // Initialize with the closest vertex to SurfacePoint (at least the first vertex will always be selected...)
    Vertex start_vertex = gc_position.nearestVertex();
    priority_queue.push({0.0, start_vertex});
    distance[start_vertex] = 0.0;

    Vector3 global_pos = getGlobalPosition();

    int iter_counter = 0;
    while (!priority_queue.empty() && iter_counter < MAX_NUMBER_OF_MESH_SEARCH_ITERS) {
        auto [dist, v] = priority_queue.top();
        priority_queue.pop();

        // If the distance exceeds the radius, stop processing
        if (dist > coverage_radius) continue;

        // Add adjacent faces of the vertex
        for (Face f : v.adjacentFaces()) {
            faces_within_radius.insert(f);
        }

        // Traverse neighbors
        for (Halfedge he : v.outgoingHalfedges()) {
            Vertex neighbor = he.tipVertex();
            //double edge_length = space->gc_geometry->edgeLengths[he.edge()];
            // TODO: decide whether it's better to use euclidean distance or the graph length (dist + edge_length)
            //  here I chose euclidean distance because the domain is flat...
            double new_dist = (global_pos - space->gc_geometry->inputVertexPositions[neighbor]).norm();

            if (distance.find(neighbor) == distance.end() || new_dist < distance[neighbor]) {
                distance[neighbor] = new_dist;
                priority_queue.push({new_dist, neighbor});
            }
        }
        iter_counter++;
    }

    // save to memory:
    if (coverage_radius == radius) {
        updated_covered_faces = true;
        memory_covered_faces = faces_within_radius;
    }
    return faces_within_radius;
}

void Agent::updateLigandReceptors(double bound_receptors, geometrycentral::Vector3 average_gradient) {
    newly_bound_receptors = bound_receptors;
    newly_average_gradient = average_gradient;
}

void Agent::stepLigandReceptors(double dt) {
    // cumulative gradient:
    double avg_gradient_norm = newly_average_gradient.norm();

    if (avg_gradient_norm > 0) {
        double c_diff = (8*radius/(3*M_PI))*avg_gradient_norm;
        double diff_LR = k_binding*c_diff*R_concentration*dt/2;

        cumulative_gradient += (diff_LR/avg_gradient_norm)*newly_average_gradient;
    }

    // Ensure stability:
    if (dt*newly_bound_receptors < R_concentration)newly_bound_receptors = R_concentration/dt;

    double d_R_concentration = k_recycled*R_star_concentration - newly_bound_receptors;
    double d_LR_concentration = newly_bound_receptors - k_internalized*LR_concentration;
    double d_R_star_concentration = k_internalized*LR_concentration - k_recycled*R_star_concentration;

    R_concentration += dt*d_R_concentration;
    LR_concentration += dt*d_LR_concentration;
    R_star_concentration += dt*d_R_star_concentration;

}

void Agent::computeNewBPRWVelocity() {
    // define Unif([0,1])
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    // Like in Pollmächer's paper:
    double gradient_following_prob = cumulative_gradient.norm()*sensitivity_to_gradient;

    double p_follow_gradient = uniform_dist(space->gen);

    if (p_follow_gradient < gradient_following_prob) {
        Vector2 new_gradient_direction;
        space->convertGlobalVectorToLocalVector(gc_position, cumulative_gradient, new_gradient_direction);

        double norm_new_dir = new_gradient_direction.norm();
        if (norm_new_dir > 0) {
            gc_local_velocity_direction = new_gradient_direction/norm_new_dir;

            // reset the cumulative gradient:
            cumulative_gradient = {0., 0., 0.};
            return;
        }
        //else you exit and rotate to random angle:
    }

    // reset cumulative gradient:
    cumulative_gradient = {0., 0., 0.};

    // select a random angle:
    double angle = (2 * M_PI)*uniform_dist(space->gen);
    gc_local_velocity_direction = gc_local_velocity_direction.rotate(angle);
}

void Agent::setLocalVelocity(std::vector<double> velocity) {
    double norm = 0.0;
    for (int d =0;d<velocity.size();d++) {
        norm += velocity[d]*velocity[d];
        gc_local_velocity_direction[d] = velocity[d];
    }
    norm = sqrt(norm);

    // update the private variables:
    speed = norm;
    gc_local_velocity_direction = gc_local_velocity_direction/norm;
}

geometrycentral::Vector2 Agent::getLocalVelocity() {
    return gc_local_velocity_direction*speed;
}
geometrycentral::Vector3 Agent::getGlobalVelocity() {
    Vector2 local_velocity = getLocalVelocity();

    Vector3 global_velocity;
    space->convertLocalVectorToGlobalVector(gc_position, local_velocity, global_velocity);
    return global_velocity;
}

geometrycentral::Vector3 Agent::getGlobalPosition() {
    return gc_position.interpolate(space->gc_geometry->inputVertexPositions);
}


