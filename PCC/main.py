from para import *
from .const import *
from .Gaussian import *
from .geom import *
from .sym import *
from .mist import check_os, print_header, get_alphabet
from .util import cmd
import logging
import matplotlib.pyplot as plt

class ePCC():

    def __init__(self):

        print_header()
        check_os()

        return

    def __read_stru(self, path, sg):
        '''
        load whole strucutre from `path`.
        :param path: structure file path
        :return:
        '''
        geom = geomParser(path)
        geom.read_geom()
        geom.log_geom_info()

        self.Natom = geom.Natom
        self.q = geom.q
        self.ele_list = geom.ele_list
        self.lattice_const = geom.lattice_const
        self.lattice_para = geom.lattice_para
        self.val = geom.val

        # find symmetric atom
        sym = symParser(sg)
        sym.igeom(self.q, self.lattice_const)
        sym.judge_sym()
        self.sym_list = sym.sym_list

        geom.judge_mol()
        self.mol_list = geom.mol_list

        sym.judge_sym_mol(self.mol_list)
        self.sym_mol = sym.sym_mol
        self.Tmol = sym.Tmol
        logging.info(self.sym_mol)

        return

    def __init_calc(self):
        logging.info('Initial calculation... ')
        if not os.path.exists('0'):
            os.mkdir('0')

        # s = open('0/gau.sh', 'w')
        cmd('yes | cp -r d/* 0/')
        cmd('cd 0 && bash gau.sh > ePCC.log')
        # s.write('cd 0\n')
        self.chg = self.__read_chg(0)

        return

    def __write_shell(self):

        if not os.path.exists('d'):
            os.mkdir('d')

        self.sym_mol_name = [get_alphabet(i) for i in range(self.Tmol)]

        s = open('d/gau.sh', 'w')
        gau = gau_wrapper()

        for i in range(self.Tmol):
            imol = self.sym_mol[i][0]
            atom_list = self.mol_list[imol]
            q_mol = [list(self.q[:, i]) for i in atom_list]
            q_mol = np.transpose(q_mol)
            e_mol = [self.ele_list[i] for i in atom_list]
            net_charge = sum([self.val[i] for i in atom_list])
            gau_file = self.sym_mol_name[i]
            gau.write_gjf(path='d', file=gau_file, q=q_mol, ele_list=e_mol, net_charge=net_charge)
            gau.write_mwfn_RESP(path='d', ele_list=e_mol, sym_file=gau_file)

            # write symmetry atoms
            sym_file_mwfn = open('d/sym_atom_{}'.format(gau_file), 'w')
            for sym in range(len(self.sym_list)):
                sym_record = []
                k = 0
                for atom in atom_list:
                    k += 1
                    if atom in self.sym_list[sym]:
                        sym_record.append(k)
                        # print(self.sym_list[sym])
                        # print(atom_list)
                        # print(sym_record)
                        # print('')
                sym_record = list(set(sym_record))
                if len(sym_record) > 1:
                    for atom in sym_record:
                        sym_file_mwfn.write('{}, '.format(atom))
                    sym_file_mwfn.write('\n')
            sym_file_mwfn.close()

            s.write(
                '''
gaucov=`grep -c -s 'Normal termination' {}.log`
if [[ "${{gaucov}}" -eq 0 ]]; then
    g16 < {}.gjf 1> {}.log 2>&1
else
    echo Gaussian calculation already completed!
fi
                '''.format(gau_file, gau_file, gau_file)
            )
            s.write(
                '''
gaucov=`grep -c -s 'Normal termination' {}.log`
if [[ "${{gaucov}}" -eq 0 ]]; then
    echo Gaussian not converged, exit!
    exit
else
    formchk {}.chk {}.fchk
fi
                '''.format(gau_file, gau_file, gau_file)
            )
            s.write(
                '''
chgcov=`grep -c -s 'If outputting atom coordinates' {}_resp.log`
if [[ "${{chgcov}}" -eq 0 ]]; then
    Multiwfn {}.chk < mwfn_{}.resp 1> {}_resp.log 2>&1
    antechamber -i {}.log -fi gout -o {}.mol2 -fo mol2 -pf y -c resp 1> {}_antechamber.log 2>&1
else
    echo RESP already fitted!
fi

chgcov=`grep -c -s 'If outputting atom coordinates' {}_resp.log`
if [[ "${{chgcov}}" -eq 0 ]]; then
    echo RESP fitting failure, check log please! exit!
exit; fi
                '''.format(gau_file, gau_file, gau_file, gau_file, gau_file, gau_file, gau_file, gau_file)
            )
        gau.write_pdb(path='d', file='0',
                      lattice_para=self.lattice_para,
                      q=self.q, ele_list=self.ele_list)

        s.close()
        return

    def __read_chg(self, it, ext='mol2'):

        chg_all = np.empty(self.Natom)
        chg_all[:] = np.nan

        for i in range(self.Tmol): # loop types of mol
            chg = []
            if ext == 'chg':
                chg_path = os.path.join(str(it), '{}.chg'.format(self.sym_mol_name[i]))
                chg_file = open(chg_path, 'r').readlines()
                chg = [float(i.split()[-1]) for i in chg_file]
            elif ext == 'mol2':
                chg_path = os.path.join(str(it), '{}.mol2'.format(self.sym_mol_name[i]))
                chg_file = open(chg_path, 'r').readlines()

                start = 0
                for l, line in enumerate(chg_file):
                    if line[0] == '@' and 'ATOM' in line:
                        start = l + 1
                        break

                while chg_file[start][0] != '@':
                    chg.append( float(chg_file[start].split()[-1]) )
                    start += 1

            for j in range(len(self.sym_mol[i])): # loop symmetry mol
                imol = self.sym_mol[i][j]
                atom_list = self.mol_list[imol]
                for kth, k in enumerate(atom_list): # loop atoms
                    chg_all[k] = chg[kth]

        for i in range(self.Natom):
            if np.isnan(chg_all[i]):
                chg_all[i] = val_chg[self.ele_list[i]]

        chg_path = os.path.join(str(it), '{}.chg'.format(it))
        chg_file = open(chg_path, 'w')
        for i in range(self.Natom):
            chg_file.write('{:}\t{:15.10}\n'.format(self.ele_list[i], chg_all[i]))
        chg_file.close()

        return chg_all

    def __it_calc(self, it):
        logging.info('Iteration : {}'.format(it))

        if not os.path.exists(str(it)):
            os.mkdir(str(it))

        cmd('yes | cp -r d/* {}'.format(it))

        geom = geomParser()
        geom.igeom(self.Natom, self.ele_list, self.q, self.lattice_const)
        q_super, ele_super, chg_super = geom.supercell(self.q, Nim=Nim, chg=self.chg)

        gau = gau_wrapper()

        for i in range(self.Tmol):
            imol = self.sym_mol[i][0]
            atom_list = self.mol_list[imol]
            q_mol = [list(self.q[:, i]) for i in atom_list]
            q_mol = np.transpose(q_mol)
            e_mol = [self.ele_list[i] for i in atom_list]
            net_charge = sum([self.val[i] for i in atom_list])
            gau_file = self.sym_mol_name[i]

            charge_q, charge_ele, charge = geom.exclude_q(q=q_mol, q_super=q_super, chg=chg_super, ele_list=ele_super)
            gau.write_gjf(path=str(it), file=gau_file, oldchk='../{}/{}'.format(it-1, gau_file),
                          q=q_mol, ele_list=e_mol, net_charge=net_charge,
                          charge=charge, charge_q=charge_q)
            gau.write_pdb(path=str(it), file=gau_file,
                          lattice_para=self.lattice_para,
                          q=q_mol, ele_list=e_mol, #chg=self.chg,
                          q_super=charge_q, ele_super=charge_ele, chg_super=charge)

        cmd('cd {} && bash gau.sh > ePCC.log'.format(it))

        return

    def __cov_check(self, it):

        chg0 = self.__read_chg(it-1)
        chg1 = self.__read_chg(it)
        logging.info('Sum of charge : {:15.8f}'.format(sum(chg1)))

        max_diff = np.max(abs(chg0 - chg1))
        logging.info('Max difference : {:15.8f}'.format(max_diff))
        if max_diff < max_diff_crt:
            return True
        else:
            return False

    def __plot_cov(self, it):
        chg = np.empty([self.Natom, it+1])
        chg_dif = np.empty([self.Natom, it])
        for i in range(it+1):
            chg[:, i] = self.__read_chg(i)
            if i > 0:
                chg_dif[:, i-1] = np.abs(chg[:, i] - chg[:, i-1])

        x = list(range(it+1))
        uniq_atom = []
        for i in range(self.Tmol): # loop symmetry molecule list
            imol = self.sym_mol[i][0] # select first molecule in each sym mol group
            uniq_atom += self.mol_list[imol]

        for j in uniq_atom:
            plt.plot(x, chg[j, :])
        plt.savefig('charge_evolution.svg', format='svg')
        plt.close()

        for j in uniq_atom:
            plt.plot(x[1:], chg_dif[j, :])
        plt.yscale('log')
        plt.savefig('charge_evolution_difference.svg', format='svg')

        return

    def main(self):

        self.__read_stru(stru_file, space_group)
        self.__write_shell()
        self.__init_calc()

        cov = False
        it = 1
        while not cov:
            self.__it_calc(it)
            self.chg = self.__read_chg(it)
            cov = self.__cov_check(it)
            it += 1

        if cov:
            self.__plot_cov(it-1)

        return