#!/home/kwalter/.conda/envs/my_gatk_4.0.0.0_kwalter/bin/python

"""haploidify.py 
   Reads in full path to input diploid single-sample VCF file and writes a haploid VCF in same location.
   Arg 1: VCF filepath. Arg 2: output haploid VCF filepath. 
"""

import pandas as pd
import vcf 
import allel
import scipy
import numpy as np
import sys
import numcodecs
import zarr
from Bio import AlignIO
from Bio import SeqIO
import sys
import argparse
import copy
import collections
import random

vcf_path= sys.argv[1]
out_path= sys.argv[2]

# define reader
vcf_reader = vcf.Reader(open(vcf_path, 'rb'))

# define writer
vcf_writer = vcf.Writer(open(out_path, 'w'), vcf_reader,
                    lineterminator='\n')

# define sample
samp = vcf_reader.samples[0]
samps = vcf_reader.samples

# loop over all records and haploidify the DeepVariant diploid output
# For lowQual calls, no AD is rseported, skip these sites, and set as noCall.

for record in vcf_reader:
  #record = next(vcf_reader)
  #print(record)
  
  # get the data fields of the tuple
  f_keys = list((record.samples[0]).data._fields)
  
  # copy record
  new_record = copy.deepcopy(record)
  
  #print(record)
  
  for sx in range(len(record.samples)):
    new_record.samples[sx].data = collections.namedtuple('CallData', f_keys)
   
    # get values
    f_vals = [record.samples[sx].data[vx] for vx in range(len(f_keys))]
    handy_dict = dict(zip(f_keys, f_vals))
  
    # get alleles called
    allele1 =  handy_dict['GT'].split('/')[0]
    allele2 =  handy_dict['GT'].split('/')[1]

    # For DeepVariant: calls that have been converted to allsites VCF, MIN_DP of homozygous ref block. 
    if 'MIN_DP' in f_keys: 
      # For DeepVariant calls that have been converted to an allsites VCF: if MIN_DP == 0, GT is nocall. 
      if handy_dict.get('MIN_DP') == 0: 
        handy_dict['GT']  = '.'
    
      # DeepVariant MIN_DP > 0, then this is a GVCF block of homozygous ref calls. GT is 0.    
      elif handy_dict.get('MIN_DP') > 0 & (allele1 == allele2) :
        handy_dict['GT']  = allele1
 
    # If a homozygous call, set to first allele call. 
    if allele1 == allele2 :
      #print(handy_dict['GT'])
      handy_dict['GT'] = allele1       
      
    # If Het call, and, if AD is called and length is not None. 
    elif ('AD' in f_keys) and  (handy_dict['AD']) : 
        #print('AD called')
        # length of AD field (includes the ref and all alleles called)
      len_AD = len(handy_dict['AD'])

      print(handy_dict)
      # access allele depths (AD field is length of the set of the ref call '0' and the two genotypes called
      AD1 = handy_dict['AD'][int(allele1)]
      AD2 = handy_dict['AD'][int(allele2)]

      # compare allele depths 
      if AD1 > AD2 :
        handy_dict['GT'] = allele1
      if AD1 < AD2 :
        handy_dict['GT'] = allele2
      if AD1 == AD2 :
        gt = random.choice([allele1,allele2])
        handy_dict['GT'] = gt

    # define new values to put in new_record    
    new_vals = [handy_dict[x] for x in f_keys]  

    # finally set CallData
    new_record.samples[sx].data = new_record.samples[sx].data._make(new_vals)
    
  # write to VCF
  vcf_writer.write_record(new_record)
  

vcf_writer.close()