#ifndef AGENT_TYPES_H
#define AGENT_TYPES_H

// NOTE-CONVENTION: good cells are positive, bad one ar negative as obvious
// best cells are good cells but faster
enum AgentType {
    AGENT_TYPE_EVIL_CELL = -1,
    AGENT_TYPE_GOOD_CELL = +1,
    AGENT_TYPE_BEST_CELL = +2
};

#endif
