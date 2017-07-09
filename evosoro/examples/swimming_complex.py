#!/usr/bin/python
"""

In this example we plug a different version of VoxCad in (_voxcad_land_water), which implements some additional features.

We evolve soft robots in a simple fluid environment (as in swimming_basic.py). In addition to some discrete phenotypic
traits (voxel material based on a predefined palette of materials), we evolve some continuous phenotypic traits:
the stiffness distribution and the actuation phase offset, which are real values associated to each voxel. We also
evolve the global actuation frequency.


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
import math

# Appending repo's root dir in the python path to enable subsequent imports
sys.path.append(os.getcwd() + "/../..")

from evosoro.base import Sim, Env, ObjectiveDict
from evosoro.networks import CPPN
from evosoro.softbot import Genotype, Phenotype, Population
from evosoro.tools.algorithms import ParetoOptimization
from evosoro.tools.utils import count_occurrences, make_material_tree, rescaled_positive_sigmoid, inverted_sigmoid
from evosoro.tools.checkpointing import continue_from_checkpoint


VOXELYZE_VERSION = '_voxcad_land_water'  # Let's source a different version of VoxCad, supporting water as well as the evolution of other attributes (stiffness distribution, actuation phase offsets)
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


MIN_ELASTIC_MOD = 0.01e006  # when, evolving stiffness, min elastic mod
MAX_ELASTIC_MOD = 1e006  # when, evolving stiffness, max elastic mod
MAX_FREQUENCY = 4.0  # when evolving a global actuation frequency, max frequency


def frequency_func(x):
    return MAX_FREQUENCY * 2.5 / (np.mean(1/x) + 1.5)  # SAM: inverse the additional inverse in read_write_voxelyze.py

# Swimming Parameters
FLUID_ENV = 1  # if 1, floor is disabled and swimming (drag) forces are added
RHO_FLUID = 1000.0  # water density
C_DRAG = 1.5  # fluid drag associated to a triangular facet
AGGREGATE_DRAG_COEF = 0.5 * C_DRAG * RHO_FLUID # aggregate drag coefficient

RUN_DIR = "swimming_complex_data"  # Subdirectory where results are going to be generated
RUN_NAME = "SwimmingComplex"
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

        # Let's add a first CPPN to the genotype. It dictates a continuous phenotypic trait,
        # the actuation phase of each voxel with respect to a global CPG-like sinusoidal signal
        self.add_network(CPPN(output_node_names=["phase_offset", "frequency"]))

        # Let's map this CPPN output to a VXA tag named <PhaseOffset>
        self.to_phenotype_mapping.add_map(name="phase_offset", tag="<PhaseOffset>", 
                                          func=partial(rescaled_positive_sigmoid, x_min=0, x_max=2*math.pi))

        self.to_phenotype_mapping.add_map(name="frequency", tag="<TempPeriod>", env_kws={"frequency": frequency_func})  # tag actually doesn't do anything here

        # Now adding a second CPPN, with three outputs. "shape" the geometry of the robot
        # (i.e. whether a particular voxel is empty or full),
        # if full, "muscleOrTissue" dictates whether a voxel is actuated or passive. The third output, "stiffness",
        # is another continuous attribute, the stiffness (elastic modulus) of each voxel
        # (overrides elastic mod defined in the materials palette)
        self.add_network(CPPN(output_node_names=["shape", "muscleOrTissue", "stiffness"]))

        # Once remapped from [-1,1] to [MIN_ELASTIC_MOD, MAX_ELASTIC_MOD] through "func",
        # the "stiffness" output goes directly to the <Stiffness> vxa tag as a continuous property.
        # We also pass min and max elastic mod as sub-tags of <Stiffness> (will be used by VoxCad)
        self.to_phenotype_mapping.add_map(name="stiffness", tag="<Stiffness>",
                                          func=partial(rescaled_positive_sigmoid, x_min=MIN_ELASTIC_MOD, x_max=MAX_ELASTIC_MOD),
                                          params=[MIN_ELASTIC_MOD, MAX_ELASTIC_MOD],
                                          param_tags=["MinElasticMod", "MaxElasticMod"])

        # The mapping for materials depends on both "shape" and "muscleOrTissue", through the following dependencies.
        # Basically, if "shape" is false (cppn output < 0), the material with id "0" is assigned (-> empty voxel).
        # If, instead, "shape" is true (cppn output > 0), we look at the "muscleOrTissue" output to determine
        # whether the material is actuated (id = 3) or passive (id = 1). These material IDs are refer to a
        # fixed palette of materials, currently hardcoded in tools/read_write_voxelyze.py
        self.to_phenotype_mapping.add_map(name="material", tag="<Data>", func=make_material_tree,
                                          dependency_order=["shape", "muscleOrTissue"], output_type=int)

        self.to_phenotype_mapping.add_output_dependency(name="shape", dependency_name=None, requirement=None,
                                                        material_if_true=None, material_if_false="0")

        self.to_phenotype_mapping.add_output_dependency(name="muscleOrTissue", dependency_name="shape",
                                                        requirement=True, material_if_true="3", material_if_false="1")


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
my_env.add_param("fluid_environment", FLUID_ENV, "<FluidEnvironment>")
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


