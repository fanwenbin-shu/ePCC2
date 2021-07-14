import numpy as np
import logging
from .util import find_mol_graph
from fractions import Fraction

class symParser:
    '''
    symmetry operations
    '''

    def __init__(self, pg=1):
        self.pg = pg

        self.sym_delta = 1e-4

        return

    def igeom(self, q, scale):
        '''
        Get the geometry from outside. `i` means input. Scale to fractional coordinate and shift to origin.
        :param q: coordinates in [Natom, 3]
        :param scale: lattice constant [3, 3]
        :return:
        '''
        self.q = np.dot(np.linalg.inv(scale), q)
        self.q -= 0.5
        return

    def ipg(self, pg):
        self.pg = pg
        self.__read_sym_operator()
        return

    def __read_sym_operator(self):

        logging.info('Reading symmetry operator `{}`...'.format(self.pg))
        path = r'space_group/SO_{}'.format(self.pg)
        f = open(path, 'r').readlines()

        Nop = int(f[1])
        op = np.zeros([3,3,Nop])
        sf = np.zeros([3,Nop])

        for i in range(Nop):
            for j in range(3):
                line = f[5*i + 3+j].strip().split()
                assert len(line) == 3, 'Wrong operation matrix! line: `{}`'.format(line)
                # print(line)
                op[:,j,i] = [float(x) for x in line]
            line = f[5*i + 3 + 3].split()
            if len(line) == 3:
                sf[:, i] = [float(Fraction(x)) for x in line]
            elif len(line) == 0:
                sf[:, i] = 0.0
            else:
                print('Wrong operation shift! line: `{}`'.format(line))
                print(line)
                exit()

        self.Nop = Nop
        self.op = op
        self.sf = sf

        logging.info('Point group number : {}'.format(self.pg))
        logging.info('Number of symmetry operator : {}'.format(Nop))
        logging.info('')

        return

    def __do_judge_sym(self):

        logging.info('Judging symmetry atoms ...')
        Natom = self.q.shape[-1]

        q = self.q
        op_mat = self.op
        sf_mat = self.sf
        Nop = self.Nop

        judge_record = np.zeros(Natom)
        sym_list = []

        for i in range(Natom - 1):
            if judge_record[i] > -50:
                judge_record[i] = -10
            sym_atom = [i]

            for j in range(i + 1, Natom):
                if judge_record[j] > -50 and judge_record[i] > -50:  # not judged atom
                    for op in range(Nop):
                        q_old = q[:, i]
                        q_new = np.dot(q[:, j], op_mat[:, :, op])
                        q_new += sf_mat[:, op]

                        half = 0.5# - self.sym_delta
                        for c in range(3):
                            while q_old[c] < -half:
                                q_old[c] += 1
                            while q_old[c] > 0.5:
                                q_old[c] -= 1
                            while q_new[c] < -half:
                                q_new[c] += 1
                            while q_new[c] > 0.5:
                                q_new[c] -= 1

                        q_dif = q_new - q_old
                        # print(i+1, j+1, q_old, q_new, q_dif)
                        if np.dot(q_dif, q_dif) < 1e-3:
                            # print(i, j, op)
                            sym_atom.append(j)
                            judge_record[j] = -100
                            break

            if len(sym_atom) > 1:
                sym_list.append(sym_atom)

        logging.debug('Writing symmetry atom to `sym_atom` file... ')
        f = open('sym_atom', 'w')
        for mol in sym_list:
            for atom in mol:
                f.write('{}, '.format(atom + 1))
            f.write('\n')

        self.sym_list = sym_list
        logging.debug('')

        return

    def judge_sym(self):

        self.__read_sym_operator()
        self.__do_judge_sym()

        return

    def judge_sym_mol(self, mol_list):

        sym_list = self.sym_list
        Nsym = len(sym_list)
        Nmol = len(mol_list)

        sym_mol = []

        for i in range(Nmol):
            for j in range(i+1, Nmol):
                same = False
                for a in mol_list[i]:
                    for b in mol_list[j]:
                        for s in sym_list:
                            if a in s and b in s:
                                same = True
                                break
                            else:
                                same = False
                    if same:
                      break
                if same:
                    sym_mol.append([i,j])

        self.sym_mol = find_mol_graph(sym_mol)

        self.Tmol = len(self.sym_mol)
        logging.info('Types of molecules : {}'.format(self.Tmol))

        return