#include "trainer.hpp"
#include "path.hpp"

bool operator<(const Individual& A, const Individual& B)
{
	return A.x > B.x;
}


bool operator<=(const Individual& A, const Individual& B)
{
	return A.x >= B.x;
}



param_t Trainer::Train(int numgenerations, vector<Individual> population)
{
	population.emplace_back(createIndividual({1,4,16,64,1,4,16,64}));
	population.emplace_back(createIndividual({-43,31,11,58,-4,23,43,45}));
	population.emplace_back(createIndividual({1,1,1,1,1,1,1,1}));
	population.emplace_back(createIndividual({1,2,3,4,1,2,3,4}));
	population.emplace_back(createIndividual({1,0,0,0,1,0,0,0}));
	
	for (int i = 0; i < 10; ++i)
	{
		population.emplace_back(createRandomIndividual());
	}
	
	sort(population.begin(), population.end());
	
	for (int i = 0; i < numgenerations; ++i)
	{
		cout << "Starting with generation " << i << " and best is : " << population[0].A << " with val: " << population[0].x << endl;

		vector<Individual> new_population;
		
		new_population.emplace_back(population[0]);
		new_population.emplace_back(createRandomIndividual());
		
		while (new_population.size() < population.size())
		{
            int i = random_give_priority_to_primeros(0,population.size());
                auto copyA = population[i].A;
            if (probability_of_true(0.3333))
            {
                // Do genetic algorithm
                int j = random_give_priority_to_primeros(0,population.size());
                if (i == j)
                    continue;
                auto copyB = population[j].A;
                cout << "Crossing " << copyA << " and " << copyB << endl;
                for (int k = 0; k < 4; ++k)
                {
                    swap(copyB[k],copyA[k]);
                }
                
                cout << "After crossing was: " << copyA << " and " << copyB << endl;
                
                perturb(copyA,5.0);
                perturb(copyB,5.0);
                new_population.emplace_back(createIndividual(copyA));
                new_population.emplace_back(createIndividual(copyB));

            } else
            {
               // Do simple perturbation
                perturb(copyA);
                new_population.emplace_back(createIndividual(copyA));
            }
		}
		sort(new_population.begin(), new_population.end());
		population = new_population;
		cout << "Done with generation " << i << " and best is : " << population[0].A << " with val: " << population[0].x << endl;
	}

	return population[0].A;
}

Individual Trainer::createIndividual(const param_t& a) const
{
	double average = 0;
	for (auto p : m_training_set)
	{
		p->set_parameters(a);
		double val = p->dfs_search(maxnumsecondsperthingy,2).value();
		average += val;
	}
	average /= m_training_set.size();
	Individual I(a,average);
	return I;
}

Individual Trainer::createRandomIndividual() const
{
	double maxval = 100.0;
	double minval = -20.0;
	param_t a;
	for (int i = 0; i < 8; ++i) a[i] = random_real(minval,maxval);
	return createIndividual(a);
}

void perturb(param_t& params, double perturbation)
{
	for (size_t i = 0; i < 8; ++i)
	{
        double r = random_real(-perturbation,perturbation);
		params[i] += random_real(-r,r); // I want the perturbation to tend to be smaller
	}
}
	
