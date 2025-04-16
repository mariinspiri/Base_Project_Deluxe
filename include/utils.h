//
// Created by spiccioni on 2/20/25.
//
#ifndef UTILS_H
#define UTILS_H

#include "Space.h"
#include "Agent.h"
#include "Field.h"




namespace utils {
    void saveAgentsPositionToFile(Space* space,std::vector<Agent>& agents, std::string file_name);
    void saveAgentsFacesToFile(Space* space, std::vector<Agent>& agents, std::string file_name);

    void saveFieldToFile(Space* space, Field& field, std::string file_name, std::string field_name="scalar field");


}
#endif //UTILS_H
