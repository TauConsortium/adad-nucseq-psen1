#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Step 1: Mapping spatial gene expression reads to the genome and microscope images with spaceranger
## Pipeline developed by Sarah J. Eger based on:
## Space ranger tutorial: https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/tutorials/count-ff-tutorial


# Download spaceranger, cellranger & reference before running this script
# https://support.10xgenomics.com/spatial-gene-expression/software/downloads/latest
# https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation

import pandas as pd
import os
import glob
import numpy as np
from multiprocessing import Pool

def run(cmd):
    print(cmd)
    os.system(cmd)

if __name__ == "__main__":
    # Set up the directories where you have stored your data    
    FASTQ_DIR = '/home/Morgane/Visium/Human'  # Expects path to directory with <file>.fastq.gz 
    IMAGE_DIR = '/home/Morgane/Visium/Human/images' # Expects path to directory with <file>.tif
    ALIGN_DIR = '/home/eger/projects/Brain_Visium/manual_alignments' # Expects path to directory with <file>.json
    REF_DIR = '/home/eger/software/refdata-gex-GRCh38-2020-A' # Expects path to directory with reference files
    
    # Set up the directory where you want to store your output
    OUTDIR = '/home/eger/projects/Brain_Visium/spaceranger_outputs'
    
    # create outdir
    if not os.path.isdir(OUTDIR):
        os.makedirs(OUTDIR)

    # set working directory
    os.chdir(OUTDIR)

    # create empty list for commands
    cmds = [] 

    # first slide
    slide1 = {'C364_c12': ['V21M02-314', 'A1'], 'C140_c12': ['V21M02-355', 'A1'], 
                'C363_c12': ['V21M02-314', 'D1'], 'C226_c12': ['V21M02-355', 'D1']}

    # loop thru images
    for ID in slide1.keys():
        cmd = ' '.join(['spaceranger count',
                        '--id='+ID, 
                        '--transcriptome='+REF_DIR,
                        '--fastqs='+FASTQ_DIR,
                        '--sample='+ID.split('_')[0],
                        '--image='+IMAGE_DIR+'/'+ID+'.tif',
                        '--slide='+slide1[ID][0],
                        '--area='+slide1[ID][1],
                        '--loupe-alignment='+ALIGN_DIR+'/'+slide1[ID][0]+'-'+slide1[ID][1]+'.json',
                        '--localcores=20',
                        '--localmem=300'])
        cmds.append(cmd)
    
    # second slide
    slide2 = {'C363_c1': ['c363-c1', 'A1'], 'C364_c1':['c364-c1', 'B1'], 
                'C381_c1': ['c381-c1', 'C1'], 'C381_c12': ['c381-c12', 'D1']}

    # loop thru images
    for ID in slide2.keys():
        cmd = ' '.join(['spaceranger count',
                        '--id='+ID, 
                        '--transcriptome='+REF_DIR,
                        '--fastqs='+FASTQ_DIR,
                        '--sample='+slide2[ID][0],
                        '--image='+IMAGE_DIR+'/'+slide2[ID][0].replace('-', '_')+'_10x.tif',
                        '--slide=V11N01-126',
                        '--area='+slide2[ID][1],
                        '--loupe-alignment='+ALIGN_DIR+'/V11N01-126-'+slide2[ID][1]+'.json',
                        '--localcores=20',
                        '--localmem=300'])
        cmds.append(cmd)
    
    # launch commands in parallel
    pool = Pool(processes = len(cmds))
    pool.map(run, cmds)
    pool.close()
    pool.join()
    
