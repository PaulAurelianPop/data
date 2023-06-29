from time import time
from am_tdm.animal_model import AnimalModel
from am_tdm.testday_model import TestDayModel
import pandas as pd
import numpy as np




PHENO_FILE = '/app/pheno_SIM_lact123_3trait_AM.txt'
TDM_FILE = '/app/pheno_SIM_lact123_5trait_TDM.txt'
PED_FILE = '/app/ped_SIM.txt'
GENO_FILE = '/app/geno_SIM.txt'




pheno = pd.read_csv(PHENO_FILE, sep='\t')
pheno = pheno[pheno['lactation_rank'] == 1]
geno = pd.read_csv(GENO_FILE, sep='\t')
ped = pd.read_csv(PED_FILE, sep='\t')
start = time()
am = AnimalModel(pheno, 1, [(9, 'cross')], [11], ped=ped, genomic_data=geno, use_blupf90_modules=True)
print(time() - start)
print(am.variance_estimator.duration)
print(am.heritabilities)


