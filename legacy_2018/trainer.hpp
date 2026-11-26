#pragma once

#include "digraph.hpp"

const double maxnumsecondsperthingy = 0.1;

struct Individual
{
	Individual(const param_t& a, double value) : A(a), x(value) {}
	param_t A;
	double x;
};

class Trainer
{
public:
    Trainer(const std::vector<DiGraph*> training_set) : m_training_set(training_set) {}
    param_t Train(int numgenerations, std::vector<Individual> population = std::vector<Individual>());
    std::vector<DiGraph*> m_training_set;
	
	Individual createIndividual(const param_t& a) const;
	Individual createRandomIndividual() const;
};

void perturb(param_t& params, double perturbation = 15.0);
