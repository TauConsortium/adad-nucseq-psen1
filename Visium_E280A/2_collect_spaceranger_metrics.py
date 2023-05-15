#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# compare auto and manual alignment

# metrics explained here:
# https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/metrics

import pandas as pd
import os
import numpy as np
import glob

if __name__ == "__main__":
    FASTQ_DIR = '/home/Morgane/Visium/Human'
    OUTDIR = '/home/eger/projects/Brain_Visium/spaceranger_outputs'
    METDIR = '/home/eger/projects/Brain_Visium/spaceranger_metrics'

    if not os.path.isdir(METDIR):
        os.makedirs(METDIR)
    
    #### collect metrics for manual
    samples = sorted(os.listdir(OUTDIR))
    samples.remove('auto_align')

    # empty dataframe to compile metrics
    df0 = pd.DataFrame()

    for samp in samples:
        met_file = os.path.join(OUTDIR, samp, 'outs', 'metrics_summary.csv')
        df = pd.read_csv(met_file)
        df0 = pd.concat([df0, df])

    df1 = df0.set_index('Sample ID').T
    
    outfile = os.path.join(METDIR, 'Manual_align_all_samples_metrics_summary.csv')
    df1.to_csv(outfile)

    #### collect metrics for auto
    samples = sorted(os.listdir(OUTDIR+'/auto_align'))

    # empty dataframe to compile metrics
    df0 = pd.DataFrame()

    for samp in samples:
        met_file = os.path.join(OUTDIR, 'auto_align', samp, 'outs', 'metrics_summary.csv')
        df = pd.read_csv(met_file)
        df0 = pd.concat([df0, df])

    df1 = df0.set_index('Sample ID').T
    
    outfile = os.path.join(METDIR, 'Auto_align_all_samples_metrics_summary.csv')
    df1.to_csv(outfile)



    #### collect Morgane's metrics
    samples = ['C363-c1_2', 'C364-c1_2', 'C381-c1_2', 'C381-c12',
                'C140-c12_Nova', 'C226-c12_Nova', 'C363-c12_Nova', 'C364-c12_Nova']

    df0 = pd.DataFrame()

    for samp in samples:
        met_file = os.path.join(FASTQ_DIR, samp, 'outs', 'metrics_summary.csv')
        df = pd.read_csv(met_file)
        print(df.loc[0, 'Sample ID'])
        df0 = pd.concat([df0, df])

    df1 = df0.set_index('Sample ID').T
    
    outfile = os.path.join(METDIR, 'Morgane_all_samples_metrics_summary.csv')
    df1.to_csv(outfile)
