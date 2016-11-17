import random
import time
import os
import cPickle
import numpy as np
import subprocess as sub

from evaluation import evaluate_all
from selection import pareto_selection, pareto_tournament_selection
from mutation import create_new_children_through_mutation
from logging import PrintLog, initialize_folders, make_gen_directories, write_gen_stats


class Optimizer(object):
    def __init__(self, sim, env, evaluation_func=evaluate_all):
        self.sim = sim
        self.env = env
        self.evaluate = evaluation_func
        self.start_time = None

    def elapsed_time(self, units="s"):
        if self.start_time is None:
            self.start_time = time.time()
        s = time.time() - self.start_time
        if units == "s":
            return s
        elif units == "m":
            return s / 60.0
        elif units == "h":
            return s / 3600.0

    def save_checkpoint(self, directory):
        if os.path.isfile("./" + directory + "/checkpoint.pickle"):
            # copy previous checkpoint in case something happens during i/o
            sub.call("cp " + directory + "/checkpoint.pickle " + directory + "/previous_checkpoint.pickle", shell=True)

        random_state = random.getstate()
        numpy_random_state = np.random.get_state()
        data = [self, random_state, numpy_random_state]

        try:
            with open(directory + '/checkpoint.pickle', 'wb') as handle:
                cPickle.dump(data, handle, protocol=cPickle.HIGHEST_PROTOCOL)

        except EOFError:
            # try one more time
            with open(directory + '/checkpoint.pickle', 'wb') as handle:
                cPickle.dump(data, handle, protocol=cPickle.HIGHEST_PROTOCOL)

    def run(self, *args, **kwargs):
        raise NotImplementedError


class PopulationBasedOptimizer(Optimizer):
    def __init__(self, sim, env, pop, selection_func, mutation_func):
        Optimizer.__init__(self, sim, env)
        self.pop = pop
        self.select = selection_func
        self.mutate = mutation_func
        self.autosuspended = False
        self.max_gens = None
        self.directory = None
        self.name = None
        self.num_random_inds = 0

    def run(self, max_hours_runtime=29, max_gens=3000, num_random_individuals=1, directory="tests_data", name="TestRun",
            max_eval_time=60, time_to_try_again=10, checkpoint_every=100, save_vxa_every=100, save_pareto=False,
            save_nets=False, save_lineages=False, continued_from_checkpoint=False):

        if self.autosuspended:
            sub.call("rm %s/AUTOSUSPENDED" % directory, shell=True)

        self.autosuspended = False
        self.max_gens = max_gens  # can add additional gens through checkpointing

        print_log = PrintLog()
        print_log.add_timer("evaluation")
        self.start_time = print_log.timers["start"]  # sync start time with logging

        # sub.call("clear", shell=True)

        if not continued_from_checkpoint:  # generation zero
            self.directory = directory
            self.name = name
            self.num_random_inds = num_random_individuals

            initialize_folders(self.pop, self.directory, self.name, save_nets, save_lineages=save_lineages)
            make_gen_directories(self.pop, self.directory, save_vxa_every, save_nets)
            sub.call("touch {}/RUNNING".format(self.directory), shell=True)
            self.evaluate(self.sim, self.env, self.pop, print_log, save_vxa_every, self.directory, self.name,
                          max_eval_time, time_to_try_again, save_lineages)
            self.select(self.pop)  # only produces dominated_by stats, no selection happening (population not replaced)
            write_gen_stats(self.pop, self.directory, self.name, save_vxa_every, save_pareto, save_nets)

        while self.pop.gen < max_gens:

            if self.pop.gen % checkpoint_every == 0:
                print_log.message("Saving checkpoint at generation {0}".format(self.pop.gen + 1), timer_name="start")
                self.save_checkpoint(self.directory)

            if self.elapsed_time(units="h") > max_hours_runtime:
                self.autosuspended = True
                print_log.message("Autosuspending at generation {0}".format(self.pop.gen + 1), timer_name="start")
                self.save_checkpoint(self.directory)
                sub.call("touch {0}/AUTOSUSPENDED && rm {0}/RUNNING".format(self.directory), shell=True)
                break

            self.pop.gen += 1
            print_log.message("Creating folders structure for this generation")
            make_gen_directories(self.pop, self.directory, save_vxa_every, save_nets)

            # update ages
            self.pop.update_ages()

            # mutation
            print_log.message("Mutation starts")
            new_children = self.mutate(self.pop, print_log=print_log)
            print_log.message("Mutation ends: successfully generated %d new children." % (len(new_children)))

            # combine children and parents for selection
            print_log.message("Now creating new population")
            self.pop.append(new_children)
            for _ in range(self.num_random_inds):
                print_log.message("Random individual added to population")
                self.pop.add_random_individual()
            print_log.message("New population size is %d" % len(self.pop))

            # evaluate fitness
            print_log.message("Starting fitness evaluation", timer_name="start")
            print_log.reset_timer("evaluation")
            self.evaluate(self.sim, self.env, self.pop, print_log, save_vxa_every, self.directory, self.name,
                          max_eval_time, time_to_try_again, save_lineages)
            print_log.message("Fitness evaluation finished", timer_name="evaluation")  # record total eval time in log

            # perform selection by pareto fronts
            new_population = self.select(self.pop)

            # print population to stdout and save all individual data
            print_log.message("Saving statistics")
            write_gen_stats(self.pop, self.directory, self.name, save_vxa_every, save_pareto, save_nets)

            # replace population with selection
            self.pop.individuals = new_population
            print_log.message("Population size reduced to %d" % len(self.pop))

        if not self.autosuspended:  # print end of run stats
            print_log.message("Finished {0} generations".format(self.pop.gen + 1))
            print_log.message("DONE!", timer_name="start")
            sub.call("touch {0}/RUN_FINISHED && rm {0}/RUNNING".format(self.directory), shell=True)


class ParetoOptimization(PopulationBasedOptimizer):
    def __init__(self, sim, env, pop):
        PopulationBasedOptimizer.__init__(self, sim, env, pop, pareto_selection, create_new_children_through_mutation)


class ParetoTournamentOptimization(PopulationBasedOptimizer):
    def __init__(self, sim, env, pop):
        PopulationBasedOptimizer.__init__(self, sim, env, pop, pareto_tournament_selection,
                                          create_new_children_through_mutation)
