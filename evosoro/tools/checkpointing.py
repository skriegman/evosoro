import random
import os
import cPickle
import numpy as np
import subprocess as sub


def continue_from_checkpoint(directory="tests_data", additional_gens=0, max_hours_runtime=29,
                             max_eval_time=60, time_to_try_again=10, checkpoint_every=1, save_vxa_every=1,
                             save_pareto=False, save_nets=False, save_lineages=False):

    # try:
    #     if sub.check_output(["ls", directory + "/RUNNING"]):
    #         sub.call("touch {}/DUPLICATE".format(directory), shell=True)
    # except:

    if os.path.isfile("./" + directory + "/RUNNING"):
        sub.call("touch {}/DUPLICATE".format(directory), shell=True)
        print "Duplicate job submitted"

    else:
        # save previous checkpoint with datetime stamp for reference
        sub.call("cp {0}/checkpoint.pickle {0}/$(date +%F_%R)_checkpoint.pickle".format(directory), shell=True)
        sub.call("touch {}/RUNNING".format(directory), shell=True)
        sub.call("rm -f %s/voxelyzeFiles/*" % directory, shell=True)  # remove partial results

        try:
            with open(directory + '/checkpoint.pickle', 'rb') as handle:
                [optimizer, random_state, numpy_random_state] = cPickle.load(handle)

        except:
            # something went wrong writing the checkpoint : use previous checkpoint
            with open(directory + '/previous_checkpoint.pickle', 'rb') as handle:
                [optimizer, random_state, numpy_random_state] = cPickle.load(handle)

        random.setstate(random_state)
        np.random.set_state(numpy_random_state)
        max_gens = optimizer.max_gens + additional_gens
        if optimizer.pop.gen < max_gens:
            optimizer.run(continued_from_checkpoint=True, max_hours_runtime=max_hours_runtime, max_gens=max_gens,
                          max_eval_time=max_eval_time, time_to_try_again=time_to_try_again,
                          checkpoint_every=checkpoint_every, save_vxa_every=save_vxa_every,
                          save_lineages=save_lineages, save_nets=save_nets, save_pareto=save_pareto)
