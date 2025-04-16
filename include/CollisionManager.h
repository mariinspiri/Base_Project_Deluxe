#ifndef COLLISIONMANAGER_H
#define COLLISIONMANAGER_H

#include <random>
#include <geometrycentral/surface/vector_heat_method.h>

#include "../include/Space.h"
#include "../include/Agent.h"


class CollisionManager {
public:
    CollisionManager(Space* s);
    // methods:
    void checkCollisions(std::vector<Agent> agents);
    bool fixCollisions(std::vector<Agent>& agents);


private:
    Space *space;
    std::vector<int> occupation_matrix;

    std::set<std::pair<int, int>> detected_collisions;
    std::map<std::pair<int, int>, std::vector<geometrycentral::Vector2>> resolution_velocities;
    std::map<std::pair<int, int>, double> resolution_lengths;

    geometrycentral::surface::VectorHeatMethodSolver heat_solver;
    std::mt19937 rng;
    std::uniform_real_distribution<double> dist_01;
};



#endif //COLLISIONMANAGER_H
