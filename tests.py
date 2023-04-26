from time import time
from am_tdm.animal_model import AnimalModel
from am_tdm.testday_model import TestDayModel
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

PHENO_FILE = '/home/trigopower/licenta/date-Radu-AM-TDM/pheno_SIM_lact123_3trait_AM.txt'
TDM_FILE = '/home/trigopower/licenta/date-Radu-AM-TDM/pheno_SIM_lact123_5trait_TDM.txt'
PED_FILE = '/home/trigopower/licenta/date-Radu-AM-TDM/ped_SIM.txt'
GENO_FILE = '/home/trigopower/licenta/date-Radu-AM-TDM/geno_SIM.txt'


def am_pblup_lact1():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 1]
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    am = AnimalModel(pheno, 'animal_ID', [('PHYS', 'cross')], ['305d_milk'], ped=ped, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_pblup_lact2():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 2]
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], ped=ped, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_pblup_lact3():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 3]
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], ped=ped, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_pblup_all():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    pheno['lact_1'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 1)
    pheno['lact_2'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 2)
    pheno['lact_3'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 3)
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [14, 15, 16], ped=ped, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_gblup_lact1():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 1]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], genomic_data=geno, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_gblup_lact2():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 2]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], genomic_data=geno, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_gblup_lact3():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 3]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], genomic_data=geno, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_gblup_all():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    pheno['lact_1'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 1)
    pheno['lact_2'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 2)
    pheno['lact_3'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 3)
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [14, 15, 16], genomic_data=geno, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_ssgblup_lact1():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 1]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], ped=ped, genomic_data=geno, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_ssgblup_lact2():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 2]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], ped=ped, genomic_data=geno, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_ssgblup_lact3():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 3]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], ped=ped, genomic_data=geno, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_ssgblup_all():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    geno = pd.read_csv(GENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    pheno['lact_1'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 1)
    pheno['lact_2'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 2)
    pheno['lact_3'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 3)
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [14, 15, 16], ped=ped, genomic_data=geno, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_pblup_lact1_gibbs():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 1]
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], ped=ped, estimation_method='gibbs', rounds=100000, burn_in=10000,
                     sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_pblup_lact2_gibbs():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 2]
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], ped=ped, estimation_method='gibbs', rounds=100000, burn_in=10000,
                     sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_pblup_lact3_gibbs():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 3]
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], ped=ped, estimation_method='gibbs', rounds=100000, burn_in=10000,
                     sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_pblup_all_gibbs():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    pheno['lact_1'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 1)
    pheno['lact_2'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 2)
    pheno['lact_3'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 3)
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [14, 15, 16], ped=ped, estimation_method='gibbs', rounds=100000,
                     burn_in=10000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_gblup_lact1_gibbs():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 1]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], genomic_data=geno, estimation_method='gibbs', rounds=100000,
                     burn_in=10000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_gblup_lact2_gibbs():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 2]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], genomic_data=geno, estimation_method='gibbs', rounds=100000,
                     burn_in=10000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_gblup_lact3_gibbs():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 3]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], genomic_data=geno, estimation_method='gibbs', rounds=100000,
                     burn_in=10000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_gblup_all_gibbs():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    pheno['lact_1'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 1)
    pheno['lact_2'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 2)
    pheno['lact_3'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 3)
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [14, 15, 16], genomic_data=geno, estimation_method='gibbs',
                     rounds=100000, burn_in=10000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_ssgblup_lact1_gibbs():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 1]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], ped=ped, genomic_data=geno, estimation_method='gibbs',
                     rounds=100000, burn_in=10000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_ssgblup_lact2_gibbs():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 2]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], ped=ped, genomic_data=geno, estimation_method='gibbs',
                     rounds=100000, burn_in=10000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_ssgblup_lact3_gibbs():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 3]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [11], ped=ped, genomic_data=geno, estimation_method='gibbs',
                     rounds=100000, burn_in=10000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def am_ssgblup_all_gibbs():
    pheno = pd.read_csv(PHENO_FILE, sep='\t')
    geno = pd.read_csv(GENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    pheno['lact_1'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 1)
    pheno['lact_2'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 2)
    pheno['lact_3'] = pheno['305d_milk'] * (pheno['lactation_rank'] == 3)
    start = time()
    am = AnimalModel(pheno, 1, [(9, 'cross')], [14, 15, 16], ped=ped, genomic_data=geno, estimation_method='gibbs',
                     rounds=100000, burn_in=10000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(am.variance_estimator.duration)
    print(am.heritabilities)


def tdm_pblup_lact1():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 1]
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), ped=ped, use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_pblup_lact2():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 2]
    pheno['lactation_rank'] -= 1
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), ped=ped, use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_pblup_lact3():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 3]
    pheno['lactation_rank'] -= 2
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), ped=ped, use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_pblup_all():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno.drop(pheno[pheno['days_in_milk'] <= 0].index)
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), ped=ped, use_blupf90_modules=True,
                       reml_conv=1e-8)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_gblup_lact1():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 1]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), genomic_data=geno,
                       use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_gblup_lact2():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 2]
    pheno['lactation_rank'] -= 1
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), genomic_data=geno,
                       use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_gblup_lact3():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 3]
    pheno['lactation_rank'] -= 2
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), genomic_data=geno,
                       use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_ssgblup_lact1():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 1]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), ped=ped, genomic_data=geno,
                       use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_ssgblup_lact2():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 2]
    pheno['lactation_rank'] -= 1
    geno = pd.read_csv(GENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), ped=ped, genomic_data=geno,
                       use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_ssgblup_lact3():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 3]
    pheno['lactation_rank'] -= 2
    geno = pd.read_csv(GENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), ped=ped, genomic_data=geno,
                       use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_pblup_lact1_gibbs():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 1]
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), ped=ped, estimation_method='gibbs',
                       rounds=400000, burn_in=300000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_pblup_lact2_gibbs():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 2]
    pheno['lactation_rank'] -= 1
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), ped=ped, estimation_method='gibbs',
                       rounds=400000, burn_in=300000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_pblup_lact3_gibbs():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 3]
    pheno['lactation_rank'] -= 2
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), ped=ped, estimation_method='gibbs',
                       rounds=400000, burn_in=300000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_gblup_lact1_gibbs():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 1]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), genomic_data=geno,
                       estimation_method='gibbs', rounds=400000, burn_in=300000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_gblup_lact2_gibbs():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 2]
    pheno['lactation_rank'] -= 1
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), genomic_data=geno,
                       estimation_method='gibbs', rounds=400000, burn_in=300000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_gblup_lact3_gibbs():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 3]
    pheno['lactation_rank'] -= 2
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), genomic_data=geno,
                       estimation_method='gibbs', rounds=400000, burn_in=300000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_ssgblup_lact1_gibbs():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 1]
    geno = pd.read_csv(GENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), ped=ped, genomic_data=geno,
                       estimation_method='gibbs', rounds=400000, burn_in=300000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_ssgblup_lact2_gibbs():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 2]
    pheno['lactation_rank'] -= 1
    geno = pd.read_csv(GENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), ped=ped, genomic_data=geno,
                       estimation_method='gibbs', rounds=400000, burn_in=300000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_ssgblup_lact3_gibbs():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    pheno = pheno[pheno['lactation_rank'] == 3]
    pheno['lactation_rank'] -= 2
    geno = pd.read_csv(GENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), ped=ped, genomic_data=geno,
                       estimation_method='gibbs', rounds=400000, burn_in=300000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_pblup_all_gibbs():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), ped=ped, estimation_method='gibbs',
                       rounds=400000, burn_in=300000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_gblup_all_gibbs():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), genomic_data=geno,
                       estimation_method='gibbs', rounds=400000, burn_in=300000, sampling=10, use_blupf90_modules=True)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_ssgblup_all_gibbs():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    geno = pd.read_csv(GENO_FILE, sep='\t')
    ped = pd.read_csv(PED_FILE, sep='\t')
    start = time()
    tdm = TestDayModel(pheno, 'animal_ID', 'lactation_rank', 'days_in_milk', [('PHYS', 'cross')], ['MY'],
                       dim_range=(5, 310), ped=ped, genomic_data=geno, estimation_method='gibbs', rounds=400000,
                       burn_in=300000, sampling=10)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)


def tdm_gblup_all():
    pheno = pd.read_csv(TDM_FILE, sep='\t')
    geno = pd.read_csv(GENO_FILE, sep='\t')
    geno = geno[geno['animal_ID'].isin(pheno['animal_ID'])]
    pheno = pheno[pheno['animal_ID'].isin(geno['animal_ID'])]
    start = time()
    tdm = TestDayModel(pheno, 1, 4, 5, [(6, 'cross')], [10], dim_range=(5, 310), genomic_data=geno,
                       use_blupf90_modules=True, reml_maxrounds=100000)
    print(time() - start)
    print(tdm.variance_estimator.duration)
    print(tdm.avg_var_G)
    print(tdm.avg_var_P)
    print(tdm.var_R)
    print(tdm.avg_heritabilities)

if __name__ == '__main__':
    tdm_ssgblup_all_gibbs()
