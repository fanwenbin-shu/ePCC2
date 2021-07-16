import os
import sys
from para import *
from .const import ename

class gau_wrapper():

    def __init__(self):
        p3elements = ['H', 'He',
                      'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                      'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar']
        self.p3elements = [i.upper() for i in p3elements]
        return

    def __get_sm(self, ele_list, charge):

        total = 0
        for ele in ele_list:
            total += ename.index(ele) + 1

        total += charge
        if total%2 == 0:
            sm = 1
        else:
            sm = 2

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
        f.write('{} {}\n'.format(int(net_charge), int(self.__get_sm(ele_list, net_charge)) ))

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

    def write_pdb(self, path, file, q, ele_list, chg=None, lattice_para=None, q_super=None, ele_super=None, chg_super=None):

        f = open(os.path.join(path, '{}.pdb'.format(file)), 'w')

        f.write('TITLE      ePCC2 - calculating Polarized Crystal Charge (PCC)\n')
        f.write('TITLE      Github : https://github.com/fanwenbin-shu/ePCC2\n')
        f.write('REMARK     Program author : Wenbin FAN (fanwenbin@shu.edu.cn)\n')
        if lattice_para is not None:
            assert len(lattice_para) >= 6, 'Please provide lattice parameters `a,b,c,alpha,beta,gamma`! '
            f.write('CRYST1{:9.3f}{:9.3f}{:9.3f}{:7.2f}{:7.2f}{:7.2f}\n'.format(*lattice_para[:6]))

        Natom = list(q.shape)[-1]
        assert len(ele_list) >= Natom, 'Shape does not match. {} {}'.format(q.shape, len(ele_list))
        if chg is None:
            chg = [0] * Natom
        for atom in range(Natom):
            f.write('{:6s}{:5d} {:^4s}'.format('ATOM', atom+1, ele_list[atom]))
            f.write('{:1s}{:3s} {:1s}{:4d}{:1s}   '.format('', 'MOL', '', 1, ''))
            f.write('{:8.3f}{:8.3f}{:8.3f}'.format(*q[:,atom]))
            f.write('{:6.2f}{:6.2f}          '.format(1.0, chg[atom]))
            f.write('{:>2s}{:2s}\n'.format(ele_list[atom], ''))
        f.write('TER\n')

        if q_super is not None:
            Nsuper = list(q_super.shape)[-1]
            for atom in range(Nsuper):
                f.write('{:6s}{:5d} {:^4s}'.format('ATOM', atom, ele_super[atom]))
                f.write('{:1s}{:3s} {:1s}{:4d}{:1s}   '.format('', 'CHG', '', 1, ''))
                f.write('{:8.3f}{:8.3f}{:8.3f}'.format(*q_super[:,atom]))
                f.write('{:6.2f}{:6.2f} {:>4d}     '.format(1.0, chg_super[atom], atom + 1))
                f.write('{:>2s}{:2s}\n'.format('X', ''))
            f.write('TER\n')

        f.write('END')
        f.close()

        return