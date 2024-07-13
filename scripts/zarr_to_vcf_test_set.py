# proxies for internet access on the cluster
%env http_proxy=http://proxy-default:3128
%env https_proxy=http://proxy-default:3128

import numpy as np
import pandas as pd

import gzip
import h5py
import zarr
import allel; print('scikit-allel', allel.__version__)
import ipytree

%matplotlib inline
from matplotlib_inline.backend_inline import set_matplotlib_formats
set_matplotlib_formats('retina', 'png')
import matplotlib.pyplot as plt
import seaborn as sns
sns.set() # sets seaborn default "prettyness:
sns.set_style("ticks")
# scale plots
sns.set_context("paper")
import matplotlib as mpl
scale = 0.8
d = dict([(k, v*scale) for (k, v) in sns.plotting_context('paper').items()])
d['figure.figsize'] = [5.4, 3.5]
mpl.rcParams.update(d)

# # citing and reflist macros
# %macro -q _cite https://gist.githubusercontent.com/kaspermunch/c4a64203651b054e18bae16792b828c6/raw/
# %macro -q _incite https://gist.githubusercontent.com/kaspermunch/77f78d8d88a816a74b5415ac0640d2b7/raw/
# %macro -q _reflist https://gist.githubusercontent.com/kaspermunch/977e903f6675b388998b1c3d367f5d55/raw/


def write_vcf(gt_array, nr_males_as_fake_females, vcf_file_name, male_id_file_name):

    header = """##fileformat=VCFv4.0
    ##phasing=partial
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##INFO=<ID=NA,Number=1,Type=Integer,Description="Not applicable">
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"""
    
    ref = callset['variants/REF'][loc_region]
    alt = callset['variants/ALT'][loc_region]
    qual = callset['variants/QUAL'][loc_region]
    
    filter_pass = callset['variants/FILTER_PASS'][loc_region]
    assert filter_pass.sum() == len(filter_pass)
    
    
        
    with gzip.open(vcf_file_name, 'wt') as f, open(male_id_file_name, 'w') as f2:
        print(header, file=f)
        row = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    
        for i in range(nr_males_as_fake_females//2):
            row.append(f'Female{i}')
        for i in range(gt_array.shape[1] - nr_males_as_fake_females):
            row.append(f'Male{i}')
            print(f'Male{i}', file=f2)
        print('\t'.join(row), file=f)
        
        for i in range(gt_array.shape[0]):
            row = ['X', str(positions[i]), f'snp{i}', str(ref[i]), str(alt[i][0]), str(qual[i]), 'PASS', 'NA', 'GT']
            
            for genotype, is_phased in zip(gt_test[1], gt_test[1].is_phased):
                row.append(f"{genotype[0]}{is_phased and '|' or '/'}{genotype[1]}")
            print('\t'.join(row), file=f)

callset = zarr.open_group('/home/kmt/baboondiversity/data/PG_panu3_zarr_12_03_2021/callset.zarr/chrX', mode='r')

gt_zarr = callset['calldata/GT']

pos = allel.SortedIndex(callset['variants/POS'])

loc_region = pos.locate_range(12000000, 120000000)
positions = pos[loc_region]

gt_region = allel.GenotypeArray(gt_zarr[loc_region])

samples = callset['samples'][:]

male_meta_data = pd.read_csv('/home/kmt/baboondiversity/data/metadata/metadata_only_males_species_ordering.txt', sep=" ")
male_ids = (male_meta_data
                   .loc[(male_meta_data.Species == "cynocephalus") & \
                        (male_meta_data.Origin == "Mikumi, Tanzania"), 'PGDP_ID']
                  )

loc_males = np.where(np.isin(samples, male_ids))[0]

gt_males = gt_region.take(loc_males, axis=1)
gt_males.is_phased = np.ones(gt_males.shape[:-1])

hap_males = gt_males.haploidify_samples()

def shuffle_along_axis(a, axis):
    idx = np.random.rand(*a.shape).argsort(axis=axis)
    return np.take_along_axis(a,idx,axis=axis)

male_ratio = 0.5
nr_males_as_fake_females = int(male_ratio * gt_males.shape[1] * 3 / 2)
data = np.array(list(zip(*[hap_males[:, i:i+2].tolist() for i in range(0, nr_males_as_fake_females, 2)])))


combined_data_truth = np.concatenate((data, gt_males[:, nr_males_as_fake_females//2:].values), axis=1)
gt_truth = allel.GenotypeArray(combined_data_truth, dtype='i1')
gt_truth.is_phased = np.ones(gt_truth.shape[:-1])

write_vcf(gt_truth, nr_males_as_fake_females, 'test.vcf.gz', 'test_males.txt')


shuffled_data = shuffle_along_axis(data, 2)
combined_data_test = np.concatenate((shuffled_data, gt_males[:, nr_males_as_fake_females//2:].values), axis=1)
gt_test = allel.GenotypeArray(combined_data_test, dtype='i1')
gt_test.is_phased = np.ones(gt_test.shape[:-1])
gt_test.is_phased[:, :nr_males_as_fake_females//2] = np.zeros((gt_test.shape[0], nr_males_as_fake_females//2))

write_vcf(gt_test, nr_males_as_fake_females, 'test.vcf.gz', 'truth_males.txt')
            
