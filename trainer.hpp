#pragma once

#include "digraph.hpp"

const double maxnumsecondsperthingy = 0.1;

struct Individual
{
	Individual(const ParamType& a, double value) : A(a), x(value) {}
	ParamType A;
	double x;
};

class Trainer
{
public:
    Trainer(const vector<DiGraph*> training_set) : m_training_set(training_set) {}
    ParamType Train(int numgenerations, vector<Individual> population = vector<Individual>());
    vector<DiGraph*> m_training_set;
	
	Individual createIndividual(const ParamType& a) const;
	Individual createRandomIndividual() const;
};

void perturb(ParamType& params, double perturbation = 15.0);
