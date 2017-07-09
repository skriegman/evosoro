import numpy as np
import itertools
import re
import scipy.ndimage as ndimage


def identity(x):
    return x


def sigmoid(x):
    return 2.0 / (1.0 + np.exp(-x)) - 1.0


def positive_sigmoid(x):
    return (1 + sigmoid(x)) * 0.5


def rescaled_positive_sigmoid(x, x_min=0, x_max=1):
    return (x_max - x_min) * positive_sigmoid(x) + x_min


def inverted_sigmoid(x):
    return sigmoid(x) ** -1


def neg_abs(x):
    return -np.abs(x)


def neg_square(x):
    return -np.square(x)


def sqrt_abs(x):
    return np.sqrt(np.abs(x))


def neg_sqrt_abs(x):
    return -sqrt_abs(x)


def mean_abs(x):
    return np.mean(np.abs(x))


def std_abs(x):
    return np.std(np.abs(x))


def count_positive(x):
    return np.sum(np.greater(x, 0))


def count_negative(x):
    return np.sum(np.less(x, 0))


def proportion_equal_to(x, keys):
    return np.mean(count_occurrences(x, keys))


def normalize(x):
    x -= np.min(x)
    x /= np.max(x)
    x = np.nan_to_num(x)
    x *= 2
    x -= 1
    return x


def xml_format(tag):
    """Ensures that tag is encapsulated inside angle brackets."""
    if tag[0] != "<":
        tag = "<" + tag
    if tag[-1:] != ">":
        tag += ">"
    return tag


def natural_sort(l, reverse):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key, reverse=reverse)


def find_between(string, start, end):
    start = string.index(start) + len(start)
    end = string.index(end, start)
    return string[start:end]


def replace_text_in_file(filename, replacements_dict):
    lines = []
    with open(filename) as infile:
        for line in infile:
            for original, target in replacements_dict.iteritems():
                line = line.replace(original, target)
            lines.append(line)
    with open(filename, 'w') as outfile:
        for line in lines:
            outfile.write(line)


def dominates(ind1, ind2, attribute_name, maximize):
    """Returns True if ind1 dominates ind2 in a shared attribute."""
    if maximize:
        return getattr(ind1, attribute_name) > getattr(ind2, attribute_name)
    else:
        return getattr(ind1, attribute_name) < getattr(ind2, attribute_name)


def count_occurrences(x, keys):
    """Count the total occurrences of any keys in x."""
    if not isinstance(x, np.ndarray):
        x = np.asarray(x)
    active = np.zeros_like(x, dtype=np.bool)
    for a in keys:
        active = np.logical_or(active, x == a)
    return active.sum()


def two_muscles(output_state):
    return np.greater(output_state, 0) + 3


def continuous_material(output_state, *args, **kwargs):
    return make_one_shape_only(output_state) * output_state


def discretize_material(output_state, num_materials=4, *args, **kwargs):
    """Discretize outputs into bins, one for each material."""
    bins = np.linspace(-1, 1, num=num_materials+1)
    return make_one_shape_only(output_state) * np.digitize(output_state, bins)


def make_material_tree(this_softbot, *args, **kwargs):

    mapping = this_softbot.to_phenotype_mapping
    material = mapping["material"]

    if material["dependency_order"] is not None:
        for dependency_name in material["dependency_order"]:
            for network in this_softbot:
                if dependency_name in network.graph.nodes():
                    mapping.dependencies[dependency_name]["state"] = \
                        network.graph.node[dependency_name]["state"] > 0

    if material["dependency_order"] is not None:
        for dependency_name in reversed(material["dependency_order"]):
            if mapping.dependencies[dependency_name]["material_if_true"] is not None:
                material["state"][mapping.get_dependency(dependency_name, True)] = \
                    mapping.dependencies[dependency_name]["material_if_true"]

            if mapping.dependencies[dependency_name]["material_if_false"] is not None:
                material["state"][mapping.get_dependency(dependency_name, False)] = \
                    mapping.dependencies[dependency_name]["material_if_false"]

    return make_one_shape_only(material["state"]) * material["state"]


def make_material_tree_single_muscle_patches(this_softbot, *args, **kwargs):

    mapping = this_softbot.to_phenotype_mapping
    material = mapping["material"]

    # for name, details in mapping.items():
    #     if details["dependency_order"] is not None:
    for dependency_name in material["dependency_order"]:
        for network in this_softbot:
            if dependency_name in network.graph.nodes():
                mapping.dependencies[dependency_name]["state"] = \
                    network.graph.node[dependency_name]["state"] > 0

    # for name, details in mapping.items():
    #     if details["dependency_order"] is not None:
    for dependency_name in reversed(material["dependency_order"]):
        if mapping.dependencies[dependency_name]["material_if_true"] is not None:
            tmpState = mapping.get_dependency(dependency_name, True)
            if dependency_name == "muscleType":
                tmpState = make_one_shape_only(tmpState)
            material["state"][tmpState] = mapping.dependencies[dependency_name]["material_if_true"]

        if mapping.dependencies[dependency_name]["material_if_false"] is not None:
            tmpState = mapping.get_dependency(dependency_name, False)
            if dependency_name == "muscleType":
                tmpState = make_one_shape_only(tmpState)
                material["state"][ndimage.morphology.binary_dilation(tmpState)] = "1"
                # print "tmpState:"
                # print tmpState
                # print "dilated:"
                # print ndimage.morphology.binary_dilation(tmpState)
            material["state"][tmpState] = mapping.dependencies[dependency_name]["material_if_false"]

    # return details["state"]
    return make_one_shape_only(material["state"]) * material["state"]


def make_one_shape_only(output_state, mask=None):
    """Find the largest continuous arrangement of True elements after applying boolean mask.

    Avoids multiple disconnected softbots in simulation counted as a single individual.

    Parameters
    ----------
    output_state : numpy.ndarray
        Network output

    mask : bool mask
        Threshold function applied to output_state

    Returns
    -------
    part_of_ind : bool
        True if component of individual

    """
    if mask is None:
        def mask(u): return np.greater(u, 0)

    # print output_state
    # sys.exit(0)

    one_shape = np.zeros(output_state.shape, dtype=np.int32)

    if np.sum(mask(output_state)) < 2:
        one_shape[np.where(mask(output_state))] = 1
        return one_shape

    else:
        not_yet_checked = []
        for x in range(output_state.shape[0]):
            for y in range(output_state.shape[1]):
                for z in range(output_state.shape[2]):
                    not_yet_checked.append((x, y, z))

        largest_shape = []
        queue_to_check = []
        while len(not_yet_checked) > len(largest_shape):
            queue_to_check.append(not_yet_checked.pop(0))
            this_shape = []
            if mask(output_state[queue_to_check[0]]):
                this_shape.append(queue_to_check[0])

            while len(queue_to_check) > 0:
                this_voxel = queue_to_check.pop(0)
                x = this_voxel[0]
                y = this_voxel[1]
                z = this_voxel[2]
                for neighbor in [(x+1, y, z), (x-1, y, z), (x, y+1, z), (x, y-1, z), (x, y, z+1), (x, y, z-1)]:
                    if neighbor in not_yet_checked:
                        not_yet_checked.remove(neighbor)
                        if mask(output_state[neighbor]):
                            queue_to_check.append(neighbor)
                            this_shape.append(neighbor)

            if len(this_shape) > len(largest_shape):
                largest_shape = this_shape

        for loc in largest_shape:
            one_shape[loc] = 1

        return one_shape


def count_neighbors(output_state, mask=None):
    """Count neighbors of each 3D element after applying boolean mask.

    Parameters
    ----------
    output_state : numpy.ndarray
        Network output

    mask : bool mask
        Threshold function applied to output_state

    Returns
    -------
    num_of_neighbors : list
        Count of True elements surrounding an individual in 3D space.

    """
    if mask is None:
        def mask(u): return np.greater(u, 0)

    presence = mask(output_state)
    voxels = list(itertools.product(*[range(x) for x in output_state.shape]))
    num_neighbors = [0 for _ in voxels]

    for idx, (x, y, z) in enumerate(voxels):
        for neighbor in [(x+1, y, z), (x-1, y, z), (x, y+1, z), (x, y-1, z), (x, y, z+1), (x, y, z-1)]:
            if neighbor in voxels:
                num_neighbors[idx] += presence[neighbor]

    return num_neighbors
