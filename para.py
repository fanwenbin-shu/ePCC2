# path of structure file
stru_file = 'CONTCAR'
# number of images in each axis
Nim = [1, 1, 1]
# space group
space_group = 113

# charge, for Gaussian calculation
charge = 5
# spin multiplicity
sm = 1
# Gaussian route card
gau_memory = 4 # unit : GB
gau_nprocs = 8
gau_method = 'b3lyp'
gau_basis = '6-31g*'
gau_maxcycle = 256
gau_conv = 8 # 8 is default
gau_optwave = False # True for stable=opt
gau_symm = False

# convergence criteria of PCC
max_diff_crt = 1e-4

# atomic charge
val_chg = {
'Mn' : 2.0, 
'Br' : -1.0, 
}
