import random
import cPickle
import numpy as np
import subprocess as sub


def continue_from_checkpoint(directory="tests_data", max_hours_runtime=29.5, max_eval_time=60, time_to_try_again=10,
                             checkpoint_every=1, additional_gens=0, num_random_individuals=1, save_lineages=False):
    try:
        if sub.check_output(["ls", directory + "/RUNNING"]):
            sub.call("touch {}/DUPLICATE".format(directory), shell=True)
    except:
        # save previous checkpoint with datetime stamp for reference
        sub.call("cp {0}/checkpoint.pickle {0}/$(date +%F_%R)_checkpoint.pickle".format(directory), shell=True)
        sub.call("touch {}/RUNNING".format(directory), shell=True)
        sub.call("rm -f %s/voxelyzeFiles/*" % directory, shell=True)  # remove partial results
        with open(directory + '/checkpoint.pickle', 'rb') as handle:
            [optimizer, random_state, numpy_random_state] = cPickle.load(handle)

        random.setstate(random_state)
        np.random.set_state(numpy_random_state)
        max_gens = optimizer.max_gens + additional_gens
        if optimizer.pop.gen < max_gens:
            optimizer.run(max_hours_runtime=max_hours_runtime, max_gens=max_gens, directory=directory,
                          max_eval_time=max_eval_time, time_to_try_again=time_to_try_again,
                          checkpoint_every=checkpoint_every, continued_from_checkpoint=True,
                          save_lineages=save_lineages, num_random_individuals=num_random_individuals)
