#!/usr/bin/python

import random
import subprocess as sub
import numpy as np
import os
import sys

# Appending repo's root dir in the python path to enable subsequent imports
sys.path.append(os.getcwd()+"/../..")

from evosoro.base import Sim, Env, ObjectiveDict
from evosoro.networks import CPPN
from evosoro.softbot import Genotype, Phenotype, Population
from evosoro.tools.algorithms import ParetoOptimization
from evosoro.tools.utils import positive_sigmoid, mean_abs, std_abs, count_negative, count_positive


sub.call("cp ../_voxcad/voxelyzeMain/voxelyze .", shell=True)
# sub.call("cp ../_voxcad/qhull .", shell=True)
# sub.call("chmod 755 ./qhull", shell=True)  # Execution right for qhull

NUM_RANDOM_INDS = 1
MAX_GENS = 1000
POPSIZE = 10
IND_SIZE = (5, 5, 4)
SIM_TIME = 10
INIT_TIME = 0.5
DT_FRAC = 0.5
GROWTH_AMPLITUDE = 0.6
MIN_TEMP_FACT = 0.2
SAVE_VXA_EVERY = 10
TIME_TO_TRY_AGAIN = 10
MAX_EVAL_TIME = 60

SEED = 1
random.seed(SEED)
np.random.seed(SEED)


class MyGenotype(Genotype):
    def __init__(self):
        Genotype.__init__(self, orig_size_xyz=IND_SIZE)
        self.add_network(CPPN(output_node_names=["initial_size"]))
        self.to_phenotype_mapping.add_map(name="initial_size", tag="<InitialVoxelSize>",
                                          logging_stats=[mean_abs, std_abs, count_negative, count_positive])

        self.add_network(CPPN(output_node_names=["start_growth_time"]))
        self.to_phenotype_mapping.add_map(name="start_growth_time", tag="<StartGrowthTime>", func=positive_sigmoid,
                                          logging_stats=[mean_abs, std_abs])

        self.add_network(CPPN(output_node_names=["growth_time"]))
        self.to_phenotype_mapping.add_map(name="growth_time", tag="<GrowthTime>", func=positive_sigmoid,
                                          logging_stats=[mean_abs, std_abs])

        self.add_network(CPPN(output_node_names=["final_size"]))
        self.to_phenotype_mapping.add_map(name="final_size", tag="<FinalVoxelSize>",
                                          logging_stats=[mean_abs, std_abs, count_negative, count_positive])

# set simulation
my_sim = Sim(dt_frac=DT_FRAC, simulation_time=SIM_TIME, min_temp_fact=MIN_TEMP_FACT, fitness_eval_init_time=INIT_TIME)

# set environment
my_env = Env()
my_env.add_param("growth_amplitude", GROWTH_AMPLITUDE, "<GrowthAmplitude>")

# set objectives for optimization
my_objective_dict = ObjectiveDict()
my_objective_dict.add_objective(name="fitness", maximize=True, tag="<NormFinalDist>")
my_objective_dict.add_objective(name="age", maximize=False, tag=None)

# initialize a pop of SoftBots
my_pop = Population(my_objective_dict, MyGenotype, Phenotype, pop_size=POPSIZE)

# set optimization procedure
my_optimization = ParetoOptimization(my_sim, my_env, my_pop)

if __name__ == "__main__":
    my_optimization.run(max_gens=MAX_GENS, save_vxa_every=SAVE_VXA_EVERY, num_random_individuals=NUM_RANDOM_INDS,
                        time_to_try_again=TIME_TO_TRY_AGAIN, max_eval_time=MAX_EVAL_TIME)
