//
// Created by spiccioni on 3/11/25.
//
#include "../include/Field.h"

Field::Field(Space *s) : space(s),
                        field(s->finite_element_space),
                        source_term(s->finite_element_space) {}

mfem::GridFunction &Field::getField() {
    return field;
}


// TODO: implement your own PDE in weak form:
void Field::setup(double dt, double diffusion_coefficient, double decay_coefficient){
    // NOTE: if you have more than one field you can think of defining these matrices only once in the Space class
    //       I chose to put them here just for tidiness

    field = 0.0;

    mfem::BilinearForm mass_form(space->finite_element_space);
    mass_form.AddDomainIntegrator(new mfem::MassIntegrator());
    mass_form.Assemble();
    mass_form.Finalize();

    // mass_matrix (M): (\int_\Omega \phi_j\phi_k)_{jk}
    mass_matrix = mass_form.SpMat();

    mfem::BilinearForm stiffness_form(space->finite_element_space);
    stiffness_form.AddDomainIntegrator(new mfem::DiffusionIntegrator());
    stiffness_form.Assemble();
    stiffness_form.Finalize();

    // stiffness_matrix (K): (\int_\Omega \Delta\phi_j\cdot\Delta\phi_k)_{jk}
    mfem::SparseMatrix stiffness_matrix = stiffness_form.SpMat();

    // now we define our specific PDE in weak form
    // here I decided to use Backward Euler method which turns the weak form of a reaction-diffusion equation into:
    // (M + dt*(D*K + \lambda M)) c(t+dt) = M*c(t) + Sources - Sinks
    //      l.h.s.          c(t+dt) = r.h.s.
    // NOTE: in this way the l.h.s is constant if dt is kept fixed and the domain doesn't change

    // I define the reaction diffusion matrix as (M + dt*(D*K + \lambda M))
    // this matrix remains constant for the whole simulation...
    reaction_diffusion_matrix = mass_matrix;
    reaction_diffusion_matrix.Add(dt*diffusion_coefficient, stiffness_matrix);
    reaction_diffusion_matrix.Add(dt*decay_coefficient, mass_matrix);

    // set the operator: to solve the linear system we use Conjugate Gradient with Gauss-Seidel's preconditioner
    solver.SetRelTol(1e-8);
    solver.SetMaxIter(50);
    solver.SetPrintLevel(0);

    // initialize sinks:
    sink_matrix = mfem::SparseMatrix(mass_matrix.Height(), mass_matrix.Width());
}

void Field::computeSources(std::vector<Agent> &agents, double source_intensity) {
    // here we define S(x,t) as an actual function, and then we project it onto the functional space (expensive but more precise)
    // if for some reason you want also you want to make the evil agents move, consider a similar implementation to computeSinkTerm()
    class EvilSources : public mfem::Coefficient {
    private:
        std::vector<Agent> agents;
        double S_0 = 0.0;

    public:
        EvilSources(std::vector<Agent> _agents, double S_0) : agents(_agents), S_0(S_0) {}

        double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip) override {
            mfem::Vector trans;
            T.Transform(ip, trans);

            for (int s=0; s!=agents.size(); s++) {
                auto &agent = agents[s];
                if (agent.getAgentType() < 0) {
                    auto evil_position = agent.getGlobalPosition();

                    double dx = trans(0) - evil_position.x;
                    double dy = trans(1) - evil_position.y;
                    double dz = trans(2) - evil_position.z;

                    double dist = sqrt(dx*dx + dy*dy + dz*dz);
                    if (dist < agent.getAgentRadius()) {
                        return S_0;
                    }
                }
            }
            return 0.0;
        }
    };

    source_term = 0.0;

    EvilSources sources(agents, source_intensity);

    source_term.AddDomainIntegrator(new mfem::DomainLFIntegrator(sources));
    source_term.Assemble();
}

void Field::computeSinksAndBindReceptors(double dt, std::vector<Agent> &agents) {
    // Here we do not define the function analytically, but construct the sink_matrix from the mass_matrix directly
    // (much much faster, compared to projecting the function)
    // Moreover since we are iterating over the faces and calculating the absorbed chemokines, we also update the bound receptors
    // inside of the agents...
    // set to 0:
    sink_matrix = 0.0;

    // TODO: also extract the average gradient and the bound-receptors

    for (auto &agent: agents) {
        if (agent.getAgentType() < 0) {continue;}
        // only good agents:
        double radius = agent.getAgentRadius();
        geometrycentral::Vector3 agent_position = agent.getGlobalPosition();

        std::unordered_set<geometrycentral::surface::Face> covered_elements;
        // filter out non-covered faces:
        for (auto face :  agent.findFacesWithinRadius()) {
            geometrycentral::Vector3 center_of_mass(0.0, 0.0, 0.0);
            for (auto vert : face.adjacentVertices()) {
                center_of_mass += space->gc_geometry->inputVertexPositions[vert];
            }
            center_of_mass /= 3;

            if ((center_of_mass - agent_position).norm() <= radius) {
                covered_elements.insert(face);
            }
        }

        double area = 0.0; // we will calculate the "covered" area which is normally a bit bigger than the actual area...

        // variables to update the agent:
        double newly_bound_receptors = 0.0;
        geometrycentral::Vector3 average_gradient{0., 0., 0.};

        // We find the Degrees Of Freedom for each element aka indexes of interested basis functions
        std::unordered_set<int> agent_dofs;
        for (auto element: covered_elements) {
            mfem::Array<int> dofs;
            int element_idx = element.getIndex();
            space->finite_element_space->GetElementDofs(element_idx, dofs);

            for (int i =0;i != dofs.Size();i++) {
                agent_dofs.insert(dofs[i]);
            }

            // Field -> Agent interaction:
            // calculate the locally bound receptors
            double element_area = space->gc_geometry->faceArea(element);
            area += element_area;

            double ligand = space->gc_geometry->faceArea(element) * field.GetValue(element_idx, {1./3., 1./3., 1./3.});
            if (!std::isinf(ligand)) {newly_bound_receptors += ligand;}

            // calculate the gradient (in global coordinates):
            mfem::ElementTransformation* trans = space->mfem_mesh->GetElementTransformation(element_idx);
            mfem::Vector gradient;
            field.GetGradient(*trans, gradient);

            for (int d =0;d != 3;d++){average_gradient[d] += gradient[d];}
        }
        average_gradient /= covered_elements.size();

        double alpha_m = agent.k_binding*agent.getFreeReceptors()/(M_PI * radius * radius);
        newly_bound_receptors *= alpha_m;

        // save the bound receptors in the agent, with the average gradient:
        agent.updateLigandReceptors(newly_bound_receptors, average_gradient);
        agent.stepLigandReceptors(dt);

        // Now that we know which functional basis we have to use we assemble the sink matrix from the mass one:
        for (int row : agent_dofs) {
            // now we extract all the non-zero elements for each row e.g. for k-th row: {(j, <value>) | \int_\Omega \phi_k\phi_j \neq 0}
            mfem::Array<int> cols; // indexes
            mfem::Vector     vals; // values

            int number_of_non_zero_els = mass_matrix.GetRow(row, cols, vals); //returns the size

            for (int j = 0; j < number_of_non_zero_els; j++){
                int   col = cols[j];
                double val = vals[j];

                // only insert into sink_matrix if col is also in the region's union_dofs
                if (agent_dofs.contains(col))
                {
                    // Add alpha * val to sink_matrix(row, col)
                    sink_matrix.Add(row, col, alpha_m * val);
                }
            }
        }
    }

}

// IMEX SCHEME: explicit sinks and sources
void Field::step(double dt) {
    mfem::Vector rhs_vector(mass_matrix.Height());

    // assemble the r.h.s. in this case: M*c + dt*F
    mass_matrix.Mult(field, rhs_vector);
    rhs_vector.Add(dt, source_term);

    mfem::Vector sink_term(mass_matrix.Height());
    sink_matrix.Mult(field, sink_term);

    rhs_vector.Add(-dt, sink_term);
    // run Backward-Euler and overwrite field:
    solver.Mult(rhs_vector, field);
}

/*
// ALTERNATIVE IMEX scheme (implicit sinks/ explicit sources)
void Field::step(double dt) {
    // now we assemble the lhs:
    mfem::SparseMatrix lhs_matrix = reaction_diffusion_matrix;
    lhs_matrix.Add(dt, sink_matrix);

    // to solve the linear system we use Conjugate Gradient with Gauss-Seidel's preconditioner
    preconditioner.SetOperator(lhs_matrix);
    solver.SetOperator(lhs_matrix);
    solver.SetPreconditioner(preconditioner);

    mfem::Vector rhs_vector(mass_matrix.Height());

    // assemble the r.h.s. in this case: M*c + dt*F
    mass_matrix.Mult(field, rhs_vector);
    rhs_vector.Add(dt, source_term);

    // run Backward-Euler and overwrite field:
    solver.Mult(rhs_vector, field);
}
*/

