from ase.ga.data import DataConnection
from ase.ga.element_mutations import RandomElementMutation
from ase.ga.element_crossovers import OnePointElementCrossover
from ase.ga.offspring_creator import OperationSelector
from ase.ga.population import Population
from ase.ga.convergence import GenerationRepetitionConvergence

from ga_fcc_alloys_relax import relax

# Specify the number of generations this script will run
num_gens = 40

db = DataConnection('fcc_alloys.db')
ref_db = 'refs.db'

# Retrieve saved parameters
population_size = db.get_param('population_size')
metals = db.get_param('metals')

# Specify the procreation operators for the algorithm
# Try and play with the mutation operators that move to nearby
# places in the periodic table
oclist = ([1, 1], [RandomElementMutation(metals),
                   OnePointElementCrossover(metals)])
operation_selector = OperationSelector(*oclist)

# Pass parameters to the population instance
pop = Population(data_connection=db,
                 population_size=population_size)

# We form generations in this algorithm run and can therefore set
# a convergence criteria based on generations
cc = GenerationRepetitionConvergence(pop, 3)

# Relax the starting population
while db.get_number_of_unrelaxed_candidates() > 0:
    a = db.get_an_unrelaxed_candidate()
    relax(a, ref_db)
    db.add_relaxed_step(a)
pop.update()

# Run the algorithm
for _ in range(num_gens):
    if cc.converged():
        print('converged')
        break
    for i in range(population_size):
        a1, a2 = pop.get_two_candidates(with_history=False)
        op = operation_selector.get_operator()
        a3, desc = op.get_new_individual([a1, a2])

        db.add_unrelaxed_candidate(a3, description=desc)

        relax(a3, ref_db)
        db.add_relaxed_step(a3)

    pop.update()

    # Print the current population to monitor the evolution
    print(['-'.join(p.get_chemical_symbols()) for p in pop.pop])
