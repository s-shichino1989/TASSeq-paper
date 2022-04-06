#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import unicode_literals

'''
Module for convert 10X BAM file to raw gene-expression count data (non UMI compression)
require pysam 0.15.4

'''
import argparse
import pandas as pd
import time
import sys
import os
import utils
import glob
import regex as re
import _version as _v
import gzip
import pysam
import csv

def main():

    des = 'BAM conversion, ' + _v.desc
    parser = argparse.ArgumentParser(description=des, add_help=True)
    
    #require cell barcode-annotated BAM
    parser.add_argument('--bam', action='store', dest='R2_bam', required=True,
                        help='10X output BAM file')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    # output and intermediate files
    R2_BAM = args.R2_bam
    final_txt = "{}_Annotated_mapping_BAM.txt.gz".format(os.path.basename(R2_BAM).split('_mapping_R2')[0])

    BAMconvert(R2_BAM, final_txt)
    
    #end = time.strftime("%Y-%m-%d %H:%M:%S")
    #print ("Finished reshaping: {}".format(end))
    return 0

def BAMconvert(R2_bam, final_txt):

    #print ("parsing BAM file...")
    fnew = gzip.open(final_txt, mode='wt', compresslevel=1)
    fbam = pysam.AlignmentFile(R2_bam, 'rb', threads=16)
    ra1 = csv.writer(fnew, delimiter='\t', lineterminator='\n')
    for rname in fbam:
        # assigned gene name
        try:
          ref_name = rname.get_tag(tag='GN')
        except:
          ref_name = '*'
            
        # cell barcode
        try:
          cell_index = rname.get_tag(tag='CB')
        except:
          cell_index = '*'
            
        # UMI
        try:
          mol_label = rname.get_tag(tag='UB')
        except:
          mol_label = '*'
        ra1.writerow([ref_name, cell_index, mol_label])
    
    fnew.close()

    return

if __name__ == '__main__':
    sys.exit(main())
