import os
import sys
from para import *
from .const import val

class gau_wrapper():

    def __init__(self):
        p3elements = ['H', 'He',
                      'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                      'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']
        self.p3elements = [i.upper() for i in p3elements]
        return

    def __get_sm(self, charge):

        if charge%2 == 0:
            sm = 2
        else:
            sm = 1

        return sm

    def write_gjf(self, path, file, q, ele_list, net_charge, antechamber=True, charge_q=None, charge=None, oldchk=None):

        Natom = q.shape[1]

        gjf_path = os.path.join(path, '{}.gjf'.format(file))
        f = open(gjf_path, 'w')

        if oldchk is not None:
            f.write('%oldchk={}.chk\n'.format(oldchk))
        f.write('%chk={}.chk\n'.format(file))
        f.write('%nprocs={}\n'.format(gau_nprocs))
        f.write('%mem={}GB\n'.format(gau_memory))
        f.write('#p {}/genecp '.format(gau_method))
        if not gau_symm:
            f.write('nosymm ')
        if gau_optwave:
            f.write('stable=opt ')
        if charge is not None:
            f.write('charge ')
        if antechamber:
            f.write('pop=mk iop(6/33=2) iop(6/42=6) ')
        f.write('\n')
        f.write('scf(maxcycle={}, conver={}) '.format(gau_maxcycle, gau_conv))
        if oldchk is not None:
            f.write('guess=read ')
        f.write('\n\n')
        f.write('PCC calculation, written by Wenbin FAN')
        f.write('\n\n')
        f.write('{} {}\n'.format(int(net_charge), int(self.__get_sm(net_charge)) ))

        for i in range(Natom):
            f.write('{}    {:18.12f}    {:18.12f}    {:18.12f}\n'.format(ele_list[i], *q[:,i]))
        f.write('\n')

        if charge is not None:
            Ncharge = len(charge)
            for i in range(Ncharge):
                f.write('{:18.12f}    {:18.12f}    {:18.12f}    {:18.12f}\n'.format(*charge_q[:,i], charge[i]))
            f.write('\n')

        # f.write('\n')
        f.write('-C -N -O -H -Cl 0\n{}\n****\n'.format(gau_basis))
        f.write('-Mn -Br 0\nSDD\n****\n\n')
        f.write('-Mn -Br 0\nSDD\n')
        f.write('\n\n\n')

        f.close()

        return

    def write_mwfn_RESP(self, ele_list, path='.', sym_file='ZZZ'):

        f = open(os.path.join(path, 'mwfn_{}.resp'.format(sym_file)), 'w')

        f.write('7\n18\n5\n1\nsym_atom_{}\n1\n'.format(sym_file))
        elements = list(set(ele_list))
        for element in elements:
            if element.upper() not in self.p3elements:
                f.write('\n')
        f.write('y\n0\n0\nq\n')
        f.close()

        return