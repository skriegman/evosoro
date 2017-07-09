#!/usr/bin/python
"""

In this example robots are allowed to develop as well as evolve. Robots can add and remove matter during their
lifetimes through linear volumetric changes.


Additional References
---------------------

This developmental model was introduced in:

    Kriegman, S., Cheney, N., Corucci, F., Bongard, J. (2017).
    A Minimal Developmental Model Can Increase Evolvability in Soft Robots.
    In Proceedings of GECCO '17, Berlin, German, July 15-19, 2017.

    Related video: https://youtu.be/gXf2Chu4L9A

"""
import random
import subprocess as sub
import numpy as np
import os
import sys

# Appending repo's root dir in the python path to enable subsequent imports
sys.path.append(os.getcwd()+"/../..")

from evosoro.base import Sim, Env, ObjectiveDict
from evosoro.networks import CPPN, DirectEncoding
from evosoro.softbot import Genotype, Phenotype, Population
from evosoro.tools.algorithms import ParetoOptimization
from evosoro.tools.utils import positive_sigmoid, mean_abs, std_abs, count_negative, count_positive
from evosoro.tools.checkpointing import continue_from_checkpoint


VOXELYZE_VERSION = '_voxcad'
# sub.call("rm ./voxelyze", shell=True)
sub.call("cp ../" + VOXELYZE_VERSION + "/voxelyzeMain/voxelyze .", shell=True)
# sub.call("chmod 755 ./voxelyze", shell=True)

NUM_RANDOM_INDS = 1
MAX_GENS = 1000
POPSIZE = 15
IND_SIZE = (5, 5, 4)
SIM_TIME = 10
INIT_TIME = 0.5
DT_FRAC = 0.5
GROWTH_AMPLITUDE = 0.5
MIN_TEMP_FACT = 0.4
SAVE_POPULATION_EVERY = 10
TIME_TO_TRY_AGAIN = 10
MAX_EVAL_TIME = 60
MAX_TIME = 0.5
SAVE_LINEAGES = False
CHECKPOINT_EVERY = 1
EXTRA_GENS = 0
RUN_DIR = "growth_data"
RUN_NAME = "Growth"


SEED = 1
random.seed(SEED)
np.random.seed(SEED)


class MyGenotype(Genotype):
    def __init__(self):
        Genotype.__init__(self, orig_size_xyz=IND_SIZE)
        self.add_network(CPPN(output_node_names=["initial_size"]))
        self.to_phenotype_mapping.add_map(name="initial_size", tag="<InitialVoxelSize>",
                                          logging_stats=[np.median, np.mean, np.std, count_negative, count_positive])

        self.add_network(CPPN(output_node_names=["final_size"]))
        self.to_phenotype_mapping.add_map(name="final_size", tag="<FinalVoxelSize>",
                                          logging_stats=[np.median, np.mean, np.std, count_negative, count_positive])

        # self.add_network(CPPN(output_node_names=["start_growth_time"]))
        # self.to_phenotype_mapping.add_map(name="start_growth_time", tag="<StartGrowthTime>", func=positive_sigmoid,
        #                                   logging_stats=[np.median, np.mean, mean_abs, np.std, std_abs])

        # self.add_network(CPPN(output_node_names=["growth_time"]))
        # self.to_phenotype_mapping.add_map(name="growth_time", tag="<GrowthTime>", func=positive_sigmoid,
        #                                   logging_stats=[np.median, np.mean, mean_abs, np.std, std_abs])

        # # The cited paper used a simpler direct encoding (no cppns):
        # self.add_network(DirectEncoding(output_node_name="initial_size", orig_size_xyz=IND_SIZE, p=0.5, scale=1))
        # self.to_phenotype_mapping.add_map(name="initial_size", tag="<InitialVoxelSize>")
        #
        # self.add_network(DirectEncoding(output_node_name="final_size", orig_size_xyz=IND_SIZE, p=0.5, scale=1))
        # self.to_phenotype_mapping.add_map(name="final_size", tag="<FinalVoxelSize>")


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
    # my_optimization.run(max_hours_runtime=MAX_TIME, max_gens=MAX_GENS, num_random_individuals=NUM_RANDOM_INDS,
    #                     directory=RUN_DIR, name=RUN_NAME, max_eval_time=MAX_EVAL_TIME,
    #                     time_to_try_again=TIME_TO_TRY_AGAIN, checkpoint_every=CHECKPOINT_EVERY,
    #                     save_vxa_every=SAVE_POPULATION_EVERY, save_lineages=SAVE_LINEAGES)

    # Here is how to use the checkpointing mechanism
    if not os.path.isfile("./" + RUN_DIR + "/pickledPops/Gen_0.pickle"):
        # start optimization
        my_optimization.run(max_hours_runtime=MAX_TIME, max_gens=MAX_GENS, num_random_individuals=NUM_RANDOM_INDS,
                            directory=RUN_DIR, name=RUN_NAME, max_eval_time=MAX_EVAL_TIME,
                            time_to_try_again=TIME_TO_TRY_AGAIN, checkpoint_every=CHECKPOINT_EVERY,
                            save_vxa_every=SAVE_POPULATION_EVERY, save_lineages=SAVE_LINEAGES)

    else:
        continue_from_checkpoint(directory=RUN_DIR, additional_gens=EXTRA_GENS, max_hours_runtime=MAX_TIME,
                                 max_eval_time=MAX_EVAL_TIME, time_to_try_again=TIME_TO_TRY_AGAIN,
                                 checkpoint_every=CHECKPOINT_EVERY, save_vxa_every=SAVE_POPULATION_EVERY,
                                 save_lineages=SAVE_LINEAGES)