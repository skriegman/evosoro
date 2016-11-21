import operator
import numpy as np
from copy import deepcopy
import math

from evosoro.networks import Network
from evosoro.tools.utils import sigmoid, xml_format, dominates


class Genotype(object):
    """A container for multiple networks, 'genetic code' copied with modification to produce offspring."""

    def __init__(self, orig_size_xyz=(6, 6, 6)):

        """
        Parameters
        ----------
        orig_size_xyz : 3-tuple (x, y, z)
            Defines the original 3 dimensions for the cube of voxels corresponding to possible networks outputs. The
            maximum number of SofBot voxel components is x*y*z, a full cube.

        """
        self.networks = []
        self.all_networks_outputs = []
        self.to_phenotype_mapping = GenotypeToPhenotypeMap()
        self.orig_size_xyz = orig_size_xyz

    def __iter__(self):
        """Iterate over the networks. Use the expression 'for n in network'."""
        return iter(self.networks)

    def __len__(self):
        """Return the number of networks in the genotype. Use the expression 'len(network)'."""
        return len(self.networks)

    def __getitem__(self, n):
        """Return network n.  Use the expression 'network[n]'."""
        return self.networks[n]

    def __deepcopy__(self, memo):
        """Override deepcopy to apply to class level attributes"""
        cls = self.__class__
        new = cls.__new__(cls)
        new.__dict__.update(deepcopy(self.__dict__, memo))
        return new

    def add_network(self, network):
        """Append a new network to this list of networks."""
        assert isinstance(network, Network)
        self.networks += [network]
        self.all_networks_outputs.extend(network.output_node_names)

    def express(self):
        """Calculate the genome networks outputs, the physical properties of each voxel for simulation"""

        for network in self:

            for name in network.graph.nodes():
                network.graph.node[name]["evaluated"] = False  # flag all nodes as unevaluated

            network.set_input_node_states(self.orig_size_xyz)  # reset the inputs

            for name in network.output_node_names:
                network.graph.node[name]["state"] = np.zeros(self.orig_size_xyz)  # clear old outputs
                network.graph.node[name]["state"] = self.calc_node_state(network, name)  # calculate new outputs

        for network in self:
            for name in network.output_node_names:
                if name in self.to_phenotype_mapping:
                    self.to_phenotype_mapping[name]["state"] = network.graph.node[name]["state"]

        for name, details in self.to_phenotype_mapping.items():
            # details["old_state"] = copy.deepcopy(details["state"])
            # SAM: moved this to mutation.py prior to mutation attempts loop
            if name not in self.all_networks_outputs:
                details["state"] = np.ones(self.orig_size_xyz, dtype=details["output_type"]) * -999
                if details["dependency_order"] is not None:
                    for dependency_name in details["dependency_order"]:
                        self.to_phenotype_mapping.dependencies[dependency_name]["state"] = None

        # create material matrix
        # for name, details in self.to_phenotype_mapping.items():
        #     details["old_state"] = copy.deepcopy(details["state"])
        #     for network in self:
        #         if name in network.output_node_names:
        #             details["state"] = network.graph.node[name]["state"]
        #         else:
        #             details["state"] = np.ones(self.orig_size_xyz, dtype=details["output_type"]) * -999
        #             if details["dependency_order"] is not None:
        #                 for dependency_name in details["dependency_order"]:
        #                     self.to_phenotype_mapping.dependencies[dependency_name]["state"] = None

        for name, details in self.to_phenotype_mapping.items():
            if details["dependency_order"] is not None:
                details["state"] = details["func"](self)

    def calc_node_state(self, network, node_name):
        """Propagate input values through the network"""
        if network.graph.node[node_name]["evaluated"]:
            return network.graph.node[node_name]["state"]

        network.graph.node[node_name]["evaluated"] = True
        input_edges = network.graph.in_edges(nbunch=[node_name])
        new_state = np.zeros(self.orig_size_xyz)

        for edge in input_edges:
            node1, node2 = edge
            new_state += self.calc_node_state(network, node1) * network.graph.edge[node1][node2]["weight"]

        network.graph.node[node_name]["state"] = new_state

        if node_name in self.to_phenotype_mapping:
            if self.to_phenotype_mapping[node_name]["dependency_order"] is None:
                return self.to_phenotype_mapping[node_name]["func"](new_state)

        return network.graph.node[node_name]["function"](new_state)


class GenotypeToPhenotypeMap(object):
    """A mapping of the relationship from genotype (networks) to phenotype (VoxCad simulation)."""

    def __init__(self):
        self.mapping = dict()
        self.dependencies = dict()

    def items(self):
        """to_phenotype_mapping.items() -> list of (key, value) pairs in mapping"""
        return [(key, self.mapping[key]) for key in self.mapping]

    def __contains__(self, key):
        """Return True if key is a key str in the mapping, False otherwise. Use the expression 'key in mapping'."""
        try:
            return key in self.mapping
        except TypeError:
            return False

    def __len__(self):
        """Return the number of mappings. Use the expression 'len(mapping)'."""
        return len(self.mapping)

    def __getitem__(self, key):
        """Return mapping for node with name 'key'.  Use the expression 'mapping[key]'."""
        return self.mapping[key]

    def __deepcopy__(self, memo):
        """Override deepcopy to apply to class level attributes"""
        cls = self.__class__
        new = cls.__new__(cls)
        new.__dict__.update(deepcopy(self.__dict__, memo))
        return new

    def add_map(self, name, tag, func=sigmoid, output_type=float, dependency_order=None, params=None, param_tags=None,
                env_kws=None, logging_stats=np.mean):
        """Add an association between a genotype output and a VoxCad parameter.

        Parameters
        ----------
        name : str
            A network output node name from the genotype.

        tag : str
            The tag used in parsing the resulting output from a VoxCad simulation.
            If this is None then the attribute is calculated outside of VoxCad (in Python only).

        func : func
            Specifies relationship between attributes and xml tag.

        output_type : type
            The output type

        dependency_order : list
            Order of operations

        params : list
            Constants dictating parameters of the mapping

        param_tags : list
            Tags for any constants associated with the mapping

        env_kws : dict
            Specifies which function of the output state to use (on top of func) to set an Env attribute

        logging_stats : func or list
            One or more functions (statistics) of the output to be logged as additional column(s) in logging

        """
        if (dependency_order is not None) and not isinstance(dependency_order, list):
            dependency_order = [dependency_order]

        if params is not None:
            assert (param_tags is not None)
            if not isinstance(params, list):
                params = [params]

        if param_tags is not None:
            assert (params is not None)
            if not isinstance(param_tags, list):
                param_tags = [param_tags]
            param_tags = [xml_format(t) for t in param_tags]

        if (env_kws is not None) and not isinstance(env_kws, dict):
            env_kws = {env_kws: np.mean}

        if (logging_stats is not None) and not isinstance(logging_stats, list):
            logging_stats = [logging_stats]

        self.mapping[name] = {"tag": xml_format(tag),
                              "func": func,
                              "dependency_order": dependency_order,
                              "state": None,
                              "old_state": None,
                              "output_type": output_type,
                              "params": params,
                              "param_tags": param_tags,
                              "env_kws": env_kws,
                              "logging_stats": logging_stats}

    def add_output_dependency(self, name, dependency_name, requirement, material_if_true=None, material_if_false=None):
        """Add a dependency between two genotype outputs.

        Parameters
        ----------
        name : str
            A network output node name from the genotype.

        dependency_name : str
            Another network output node name.

        requirement : bool
            Dependency must be this

        material_if_true : int
            The material if dependency meets pre-requisite

        material_if_false : int
            The material otherwise

        """
        self.dependencies[name] = {"depends_on": dependency_name,
                                   "requirement": requirement,
                                   "material_if_true": material_if_true,
                                   "material_if_false": material_if_false,
                                   "state": None}

    def get_dependency(self, name, output_bool):
        """Checks recursively if all boolean requirements were met in dependent outputs."""
        if self.dependencies[name]["depends_on"] is not None:
            dependency = self.dependencies[name]["depends_on"]
            requirement = self.dependencies[name]["requirement"]
            return np.logical_and(self.get_dependency(dependency, True) == requirement,
                                  self.dependencies[name]["state"] == output_bool)
        else:
            return self.dependencies[name]["state"] == output_bool


class Phenotype(object):
    """Physical manifestation of the genotype - determines the physiology of an individual."""

    def __init__(self, genotype):

        """
        Parameters
        ----------
        genotype : Genotype()
            Defines particular networks (the genome).

        """
        self.genotype = genotype
        self.genotype.express()

    def __deepcopy__(self, memo):
        """Override deepcopy to apply to class level attributes"""
        cls = self.__class__
        new = cls.__new__(cls)
        new.__dict__.update(deepcopy(self.__dict__, memo))
        return new

    def is_valid(self):
        """Ensures a randomly generated phenotype is valid (checked before adding individual to a population).

        Returns
        -------
        is_valid : bool
        True if self is valid, False otherwise.

        """
        for network in self.genotype:
            for output_node_name in network.output_node_names:
                if np.isnan(network.graph.node[output_node_name]["state"]).any():
                    return False
        return True


class SoftBot(object):
    """A SoftBot is a 3D creature composed of a continuous arrangement of connected voxels with varying softness."""

    def __init__(self, max_id, objective_dict, genotype, phenotype):

        """Initialize an individual SoftBot for physical simulation within VoxCad.

        Parameters
        ----------
        max_id : the lowest id tag unused
            An index to keep track of evolutionary history.

        objective_dict : ObjectiveDict()
            Defines the objectives to optimize.

        genotype : Genotype cls
            Defines the networks (genome).

        phenotype : Phenotype cls
            The physical manifestation of the genotype which defines an individual in simulation.

        """
        self.genotype = genotype()  # initialize new random genome
        self.phenotype = phenotype(self.genotype)  # calc phenotype from genome

        self.id = max_id
        self.md5 = "none"
        self.dominated_by = []  # other individuals in the population that are superior according to evaluation
        self.pareto_level = 0
        self.selected = 0  # survived selection
        self.variation_type = "newly_generated"  # (from parent)
        self.parent_genotype = self.genotype  # default for randomly generated ind
        self.parent_id = -1
        self.age = 0

        # set the objectives as attributes of self (and parent)
        self.objective_dict = objective_dict
        for rank, details in objective_dict.items():
            if details["name"] != "age":
                setattr(self, details["name"], details["worst_value"])
            setattr(self, "parent_{}".format(details["name"]), details["worst_value"])

    def __deepcopy__(self, memo):
        """Override deepcopy to apply to class level attributes"""
        cls = self.__class__
        new = cls.__new__(cls)
        new.__dict__.update(deepcopy(self.__dict__, memo))
        return new


class Population(object):
    """A population of SoftBots."""

    def __init__(self, objective_dict, genotype, phenotype, pop_size=30):

        """Initialize a population of individual SoftBots.

        Parameters
        ----------
        objective_dict : ObjectiveDict()
            Defines the objectives to optimize.

        genotype : Genotype
            The genetic code used to create an individual and passed to offspring (with modification).

        phenotype : Phenotype
            The physical manifestation of the genotype which defines an individual in simulation.

        pop_size : int
            The target number of individuals to maintain in the population.

        """
        self.genotype = genotype
        self.phenotype = phenotype
        self.pop_size = pop_size
        self.gen = 0
        self.total_evaluations = 0
        self.already_evaluated = {}
        self.all_evaluated_individuals_ids = []
        self.objective_dict = objective_dict
        self.best_fit_so_far = objective_dict[0]["worst_value"]
        self.individuals = []
        self.max_id = 0
        self.non_dominated_size = 0

        while len(self) < pop_size:
            self.add_random_individual()

    def __iter__(self):
        """Iterate over the individuals. Use the expression 'for n in population'."""
        return iter(self.individuals)

    def __contains__(self, n):
        """Return True if n is a SoftBot in the population, False otherwise. Use the expression 'n in population'."""
        try:
            return n in self.individuals
        except TypeError:
            return False

    def __len__(self):
        """Return the number of individuals in the population. Use the expression 'len(population)'."""
        return len(self.individuals)

    def __getitem__(self, n):
        """Return individual n.  Use the expression 'population[n]'."""
        return self.individuals[n]

    def pop(self, index=None):
        """Remove and return item at index (default last)."""
        return self.individuals.pop(index)

    def append(self, individuals):
        """Append a list of new individuals to the end of the population.

        Parameters
        ----------
        individuals : list of/or SoftBot
            A list of individual SoftBots to append or a single SoftBot to append

        """
        if type(individuals) == list:
            for n in range(len(individuals)):
                if type(individuals[n]) != SoftBot:
                    raise TypeError("Non-SoftBot added to the population")
            self.individuals += individuals

        elif type(individuals) == SoftBot:
            self.individuals += [individuals]

    def sort(self, key, reverse=False):
        """Sort individuals by their attributes.

        Parameters
        ----------
        key : str
            An individual-level attribute.

        reverse : bool
            True sorts from largest to smallest (useful for maximizing an objective).
            False sorts from smallest to largest (useful for minimizing an objective).

        """
        return self.individuals.sort(reverse=reverse, key=operator.attrgetter(key))

    def add_random_individual(self):
        valid = False
        while not valid:
            ind = SoftBot(self.max_id, self.objective_dict, self.genotype, self.phenotype)
            if ind.phenotype.is_valid():
                self.individuals.append(ind)
                self.max_id += 1
                valid = True

    def update_ages(self):
        """Increment the age of each individual."""
        for ind in self:
            ind.age += 1
            ind.variation_type = "survived"

    def sort_by_objectives(self):
        """Sorts the population multiple times by each objective, from least important to most important."""
        for ind in self:
            if math.isnan(ind.fitness):
                ind.fitness = self.objective_dict[0]["worst_value"]
                print "FITNESS WAS NAN, RESETTING IT TO:", self.objective_dict[0]["worst_value"]

        self.sort(key="id", reverse=True)  # (max) promotes neutral mutation
        self.sort(key="age", reverse=False)  # (min) protects younger, undeveloped solutions

        for rank in reversed(range(len(self.objective_dict))):
            goal = self.objective_dict[rank]
            self.sort(key=goal["name"], reverse=goal["maximize"])

        self.sort(key="pareto_level", reverse=False)  # min

        # print "rank in sort_by_objectives:", [len(i.dominated_by) for i in self]
        # print "age in sort_by_objectives:", [i.age for i in self]
        # print "fitness in sort_by_objectives:", [i.fitness for i in self]

    def dominated_in_multiple_objectives(self, ind1, ind2):
        """Calculate if ind1 is dominated by ind2 according to all objectives in objective_dict.

        If ind2 is better or equal to ind1 in all objectives, and strictly better than ind1 in at least one objective.

        """
        # losses = []  # 2 dominates 1
        wins = []  # 1 dominates 2
        for rank in reversed(range(len(self.objective_dict))):
            goal = self.objective_dict[rank]
            # losses += [dominates(ind2, ind1, goal["name"], goal["maximize"])]  # ind2 dominates ind1?
            wins += [dominates(ind1, ind2, goal["name"], goal["maximize"])]  # ind1 dominates ind2?
        # return np.any(losses) and not np.any(wins)
        return not np.any(wins)

    def calc_dominance(self):
        """Determine which other individuals in the population dominate each individual."""

        self.sort(key="id", reverse=False)  # if tied on all objectives, give preference to newer individual

        # clear old calculations of dominance
        self.non_dominated_size = 0
        for ind in self:
            ind.dominated_by = []
            ind.pareto_level = 0

        for ind in self:
            for other_ind in self:
                # if (other_ind.fitness >= ind.fitness) and \
                #         self.dominated_in_multiple_objectives(ind, other_ind) and \
                #         (ind.id not in other_ind.dominated_by) and \
                #         (other_ind.id != ind.id):
                if other_ind.id != ind.id:
                    if self.dominated_in_multiple_objectives(ind, other_ind) and (ind.id not in other_ind.dominated_by):
                        ind.dominated_by += [other_ind.id]

            if ind.fitness == self.objective_dict[0]["worst_value"]:  # extra penalty for doing nothing or being invalid
                ind.dominated_by += [ind.id for _ in range(self.pop_size * 2)]

            ind.pareto_level = len(ind.dominated_by)  # update the pareto level

            # update the count of non_dominated individuals
            if ind.pareto_level == 0:
                self.non_dominated_size += 1
