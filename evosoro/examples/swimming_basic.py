#!/usr/bin/python
"""

In this example we plug a different version of VoxCad in (_voxcad_land_water), which implements some additional
features. One of them is a simple mesh-based fluid model, that allows you to observe how soft morphologies evolve in a
fluid environment.

The setup is almost identical to basic.py (same genotype and phenotype, based on a discrete materials palette), with the
addition of few parameters involved in the fluid model.

The .vxa files generated in this example will embed two additional tags (FluidEnvironment and AggregateDragCoefficient)
telling to the physics engine that we want to simulate a fluid environment (floor and gravity are disabled, drag forces
are added).


References
---------------------

If using this code or the "_voxcad_land_water" version of the simulator, please cite:

  - F. Corucci, N. Cheney, H. Lipson, C. Laschi, and J. Bongard,
    "Evolving swimming soft-bodied creatures"
    Late breaking proceedings of The Fifteenth International Conference on the Synthesis and Simulation of Living Systems (ALIFE XV), 2016
    Available at: http://sssa.bioroboticsinstitute.it/sites/default/files/user_docs/2016_ALIFE_LateBreaking_SWIMMING.pdf
    Related video: https://youtu.be/4ZqdvYrZ3ro

"""
import random
import numpy as np
import subprocess as sub
from functools import partial
import os
import sys

# Appending repo's root dir in the python path to enable subsequent imports
sys.path.append(os.getcwd() + "/../..")

from evosoro.base import Sim, Env, ObjectiveDict
from evosoro.networks import CPPN
from evosoro.softbot import Genotype, Phenotype, Population
from evosoro.tools.algorithms import ParetoOptimization
from evosoro.tools.utils import count_occurrences, make_material_tree
from evosoro.tools.checkpointing import continue_from_checkpoint


VOXELYZE_VERSION = '_voxcad_land_water'  # Let's source a different version of VoxCad, supporting water
# sub.call("rm ./voxelyze", shell=True)
sub.call("cp ../" + VOXELYZE_VERSION + "/voxelyzeMain/voxelyze .", shell=True)  # Making sure to have the most up-to-date version of the Voxelyze physics engine
# sub.call("chmod 755 ./voxelyze", shell=True)
# sub.call("cp ../_voxcad/qhull .", shell=True)  # Auxiliary qhull executable, used in some experiments to compute the convex hull of the robot
# sub.call("chmod 755 ./qhull", shell=True)  # Execution right for qhull

NUM_RANDOM_INDS = 1  # Number of random individuals to insert each generation (increases diversity, works with AFPO when minimizing "age")
MAX_GENS = 1000  # Number of generations
POPSIZE = 10  # Population size (number of individuals in the population)
IND_SIZE = (6, 6, 6)  # Bounding box dimensions (x,y,z). e.g. IND_SIZE=(6, 6, 6) -> workspace is a cube of 6x6x6 voxels
SIM_TIME = 5.0  # evaluation time (seconds), including INIT_TIME!
INIT_TIME = 0.5 # initial transient time (robot settles down before actuation and fitness computation start)
DT_FRAC = 0.9  # Fraction of the optimal integration step. The lower, the more stable (and slower) the simulation.

TIME_TO_TRY_AGAIN = 30  # (seconds) wait this long before assuming simulation crashed and resending
MAX_EVAL_TIME = 60  # (seconds) wait this long before giving up on evaluating this individual
SAVE_LINEAGES = False
MAX_TIME = 8  # (hours) how long to wait before autosuspending
EXTRA_GENS = 0  # extra gens to run when continuing from checkpoint

# Fluid parameters (will be passed to Environment)
RHO_FLUID = 1000.0  # fluid density (water)
C_DRAG = 1.5  # fluid drag associated to a triangular facet
AGGREGATE_DRAG_COEF = 0.5 * C_DRAG * RHO_FLUID  # aggregate drag coefficient

RUN_DIR = "swimming_basic_data"  # Subdirectory where results are going to be generated
RUN_NAME = "SwimmingBasic"
CHECKPOINT_EVERY = 1  # How often to save an snapshot of the execution state to later resume the algorithm
SAVE_POPULATION_EVERY = 1  # How often (every x generations) we save a snapshot of the evolving population

SEED = 1 # Random seed
random.seed(SEED)  # Initializing the random number generator for reproducibility
np.random.seed(SEED)


# Defining a custom genotype, inheriting from base class Genotype
class MyGenotype(Genotype):
    def __init__(self):

        # We instantiate a new genotype for each individual which must have the following properties
        Genotype.__init__(self, orig_size_xyz=IND_SIZE)

        # The genotype consists of a single Compositional Pattern Producing Network (CPPN),
        # with multiple inter-dependent outputs determining the material constituting each voxel
        # (e.g. two types of active voxels, actuated with a different phase, two types of passive voxels, softer and stiffer)
        # The material IDs that you will see in the phenotype mapping dependencies refer to a predefined palette of materials
        # currently hardcoded in tools/read_write_voxelyze.py:
        # (0: empty, 1: passiveSoft, 2: passiveHard, 3: active+, 4:active-),
        # but this can be changed.
        self.add_network(CPPN(output_node_names=["shape", "muscleOrTissue", "muscleType", "tissueType"]))

        self.to_phenotype_mapping.add_map(name="material", tag="<Data>", func=make_material_tree,
                                          dependency_order=["shape", "muscleOrTissue", "muscleType", "tissueType"], output_type=int)

        self.to_phenotype_mapping.add_output_dependency(name="shape", dependency_name=None, requirement=None,
                                                        material_if_true=None, material_if_false="0")

        self.to_phenotype_mapping.add_output_dependency(name="muscleOrTissue", dependency_name="shape",
                                                        requirement=True, material_if_true=None, material_if_false=None)

        self.to_phenotype_mapping.add_output_dependency(name="tissueType", dependency_name="muscleOrTissue",
                                                        requirement=False, material_if_true="1", material_if_false="2")

        self.to_phenotype_mapping.add_output_dependency(name="muscleType", dependency_name="muscleOrTissue",
                                                        requirement=True, material_if_true="3", material_if_false="4")

# Define a custom phenotype, inheriting from the Phenotype class
class MyPhenotype(Phenotype):
    def is_valid(self, min_percent_full=0.3, min_percent_muscle=0.1):
        # override super class function to redefine what constitutes a valid individuals
        for name, details in self.genotype.to_phenotype_mapping.items():
            if np.isnan(details["state"]).any():
                return False
            if name == "material":
                state = details["state"]
                # Discarding the robot if it doesn't have at least a given percentage of non-empty voxels
                if np.sum(state>0) < np.product(self.genotype.orig_size_xyz) * min_percent_full:
                    return False
                # Discarding the robot if it doesn't have at least a given percentage of muscles (materials 3 and 4)
                if count_occurrences(state, [3, 4]) < np.product(self.genotype.orig_size_xyz) * min_percent_muscle:
                    return False
        return True


# Setting up the simulation object
my_sim = Sim(dt_frac=DT_FRAC, simulation_time=SIM_TIME, fitness_eval_init_time=INIT_TIME)

# Setting up the environment object
my_env = Env(sticky_floor=0, time_between_traces=0)
# Here we tell the physics engine that we want to simulate a fluid environment
my_env.add_param("fluid_environment", 1, "<FluidEnvironment>")
my_env.add_param("aggregate_drag_coefficient", AGGREGATE_DRAG_COEF, "<AggregateDragCoefficient>")


# Now specifying the objectives for the optimization.
# Creating an objectives dictionary
my_objective_dict = ObjectiveDict()

# Adding an objective named "fitness", which we want to maximize. This information is returned by Voxelyze
# in a fitness .xml file, with a tag named "NormFinalDist"
my_objective_dict.add_objective(name="fitness", maximize=True, tag="<normAbsoluteDisplacement>")

# Add an objective to minimize the age of solutions: promotes diversity
my_objective_dict.add_objective(name="age", maximize=False, tag=None)

# Adding another objective called "num_voxels", which we want to minimize in order to minimize
# the amount of material employed to build the robot, promoting at the same time non-trivial
# morphologies.
# This information can be computed in Python (it's not returned by Voxelyze, thus tag=None),
# which is done by counting the non empty voxels (material != 0) composing the robot.
my_objective_dict.add_objective(name="num_voxels", maximize=False, tag=None,
                                node_func=np.count_nonzero, output_node_name="material")

# Adding another objective named "energy", which should be minimized.
# This information is not returned by Voxelyze (tag=None): it is instead computed in Python.
# We also specify how energy should be computed, which is done by counting the occurrences of
# active materials (materials number 3 and 4)
my_objective_dict.add_objective(name="energy", maximize=False, tag=None,
                                node_func=partial(count_occurrences, keys=[3, 4]),
                                output_node_name="material")


# Initializing a population of SoftBots
my_pop = Population(my_objective_dict, MyGenotype, MyPhenotype, pop_size=POPSIZE)

# Setting up our optimization
my_optimization = ParetoOptimization(my_sim, my_env, my_pop)

# And, finally, our main
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


