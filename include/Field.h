//
// Created by spiccioni on 3/11/25.
//
#ifndef FIELD_H
#define FIELD_H

#include <iostream>
#include <mfem.hpp>

#include "Space.h"
#include "Agent.h"



class Field {
public:
    Field(Space *s);

    mfem::GridFunction &getField();

    void setup(double dt, double diffusion_coefficient, double decay_coefficient);
    void computeSources(std::vector<Agent> &agents, double source_intensity);
    void computeSinksAndBindReceptors(double dt, std::vector<Agent> &agents);

    void step(double dt);

private:
    mfem::GridFunction field;

    // solvers:
    mfem::CGSolver solver;
    mfem::GSSmoother preconditioner;

    // PDE terms:
    mfem::SparseMatrix mass_matrix;
    mfem::SparseMatrix reaction_diffusion_matrix;
    mfem::SparseMatrix sink_matrix;
    mfem::LinearForm source_term;

    Space *space;
};
#endif //FIELD_H
