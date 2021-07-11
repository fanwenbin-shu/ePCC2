import os
import numpy as np

def angle_between(v1, v2):
    '''
    calculate angles between two vectors
    ref: https://stackoverflow.com/a/13849249/71522
    :param v1:
    :param v2:
    :return:
    '''
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def find_mol_graph(l):
    '''
    find a whole molecule according to topology.
    ref: https://stackoverflow.com/questions/4842613/merge-lists-that-share-common-elements/4842897
    :param l: a list consisting several arrays with common elements, like [[1], [1,2], [3,4], [4,5], [6,7]]
    :return: an integrate list, like [[1,2], [3,4,5], [6,7]]
    '''
    if len(l) == 0:
        return []
    out = []
    while len(l) > 0:
        first, *rest = l
        first = set(first)

        lf = -1
        while len(first) > lf:
            lf = len(first)

            rest2 = []
            for r in rest:
                if len(first.intersection(set(r))) > 0:
                    first |= set(r)
                else:
                    rest2.append(r)
            rest = rest2

        out.append(first)
        l = rest

    out_list = []
    for item in out:
        out_list.append(sorted(item))

    return out_list

def const2para(lattice_const):
    '''
    Convert to lattice constant to lattice parameter.

    :param lattice_const: (3,3) array for (xyz, abc)
    :return: lattice constant
    '''
    assert lattice_const.shape == (3, 3), 'Wrong shape of lattice constant! '

    va = lattice_const[:, 0]  # vector a
    vb = lattice_const[:, 1]
    vc = lattice_const[:, 2]

    a, b, c = [np.linalg.norm(x) for x in [va, vb, vc]]

    alpha = angle_between(vb, vc)
    beta = angle_between(va, vc)
    gamma = angle_between(va, vb)

    alpha, beta, gamma = [x * 180 / np.pi for x in [alpha, beta, gamma]]

    return [a, b, c, alpha, beta, gamma]

def get_pbc_list(im):
    '''
    get an integer list from -im to im, like [-2, -1, 0, 1, 2] when im=2.
    :param im: integer
    :return: list
    '''
    l = []
    if im == 0:
        return [0]
    else:
        for i in range(im+1):
            if i != 0:
                l.append(-i)
                l.append(i)
            else:
                l.append(0)
    l.sort()

    return l

def cmd(command, path='.'):

    os.popen('cd {}'.format(path))
    res = os.popen(command)

    return res.read()