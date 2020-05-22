#!/usr/bin/python3

# This script to be used to run liftOver on genotype data stored in
# the plink format.
# See: http://genome.sph.umich.edu/wiki/LiftOver
# Downloaded from: http://genome.sph.umich.edu/wiki/LiftMap.py
#
# Modified by Scott Ritchie:
#  - to work with user specified chain files, rather than
#    the original developer's specific chain file.
#  - to not rely the original developer's path to liftOver.
#  - to provide helpful usage documentation.
#  - to clean up the intermediary BED files to avoid confusion.
#  - to generally be slightly more PEP compliant.
# 

import sys
import os
import argparse
import gzip
import tempfile
import pandas as pd
import numpy as np
#from string import Template
import basic_fun


def die(msg):
    sys.exit(msg)

def myopen(fn):
    try:
        h = gzip.open(fn)
        ln = h.read(2) # read arbitrary bytes so check if @param fn is a gzipped file
    except:
        # cannot read in gzip format
        return open(fn)
    h.close()
    return gzip.open(fn)

def bim2bed(fin, fout):
    print("Converting PLINK BIM file to UCSC BED file...")
    fin = pd.read_csv(fin, sep='\s+', header=None, names=['chr', 'snp', 'genetic', 'pos', 'a1', 'a2']) # a2 a1 are ref alt
    fin['pos0'] = fin['pos']-1
    fin[['chr', 'pos0', 'pos', 'snp']].to_csv(fout, sep='\t', index=False, header=False)
    return True

# global var:
LIFTED_SET = set()
UNLIFTED_SET = set()
def liftBed(oldBED, newBED, funlifted, chainFile, liftOverPath):
    print("Lifting BED file...")
    cmd = [liftOverPath, oldBED, chainFile, newBED, funlifted]
    basic_fun.subprocess_popen(cmd)
    #record lifted/unliftd rs
    for ln in myopen(funlifted):
        if len(ln) == 0 or ln[0] == '#':continue
        UNLIFTED_SET.add(ln.strip().split()[-1])
    for ln in myopen(newBED):
        if len(ln) == 0 or ln[0] == '#':continue
        LIFTED_SET.add(ln.strip().split()[-1])

    return True

def extractBfile(oldBfile, BED, newBfile):
    print("Extracting lifted variants from PLINK files...")
    temp = tempfile.NamedTemporaryFile('w+t', delete=False)
    snp = temp.name
    for line in open(BED):
        line = line.strip().split()
        temp.write(line[3]+'\n')
    temp.close()
    basic_fun.subprocess_popen(['plink', '--bfile', oldBfile, '--extract', snp, '--memory', args.memory, '--make-bed', 
        '--out', newBfile])
    return True

def updatebim(bim, bed):
    print("Updating SNP info in bim to new Assembly...")
    bim = pd.read_csv(bim, sep='\s+', header=None, names=['chr', 'snp', 'genetic', 'pos', 'a1', 'a2'])
    bim = bim[['snp', 'genetic', 'a1', 'a2']]
    bed = pd.read_csv(bed, sep='\s+', header=None, names=['chrom', 'pos0', 'pos', 'snp'])
    bed['chr'] = bed['chrom'].str.replace('^chr', '') # str 
    bed = bed.merge(bim, on='snp')
    bed[['chr', 'snp', 'genetic', 'pos', 'a1', 'a2']].to_csv(bim, sep='\t', index=False, header=False)
    return True

def liftDat(fin, fout):
    fo = open(fout, 'w')
    for ln in myopen(fin):
        if len(ln) == 0 or ln[0] != 'M':
            fo.write(ln)
        else:
            t, rs = ln.strip().split()
            if rs in LIFTED_SET:
                fo.write(ln)
    fo.close()
    return True

def makesure(result, succ_msg, fail_msg = "ERROR"):
    if result:
        print('SUCC: ', succ_msg)
    else:
        print('FAIL: ', fail_msg)
        sys.exit(2)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="%(prog)s converts genotype data stored in plink's PED+MAP " +
                    "format from one genome build to another, using liftOver."
    )
    parser.add_argument('-b', "--bfile", dest='bfile', required = True,
                        help='The plink BED and BIM files to `liftOver`.')
    parser.add_argument('-o', "--out", dest='prefix', required = True,
                        help='The prefix to give to the output files.')
    parser.add_argument('-c', "--chain", dest='chainFile', required = True,
                        help='The location of the chain file to provide to ' +
                             '`liftOver`.')
    parser.add_argument('-e', "--bin", dest='liftOverExecutable', default='liftOver',
                        help='The location of the `liftOver` executable.')
    parser.add_argument('-m', '--memory', dest='memory', default='10000',
                        help='memory for PLINK')
    parser.add_argument('-d', "--dat", dest='datFile',
                        help='[NOT USED]Optionally remove "unlifted SNPs" from a data ' +
                             'file containing a list of SNPs (e.g. for ' +
                             ' --exclude or --include in `plink`)')

    # Show usage message if user hasn't provided any arguments, rather
    # than giving a non-descript error message with the usage()
    if len(sys.argv) == 1:
      parser.print_help()
      sys.exit()

    args = parser.parse_args()

    # make BED file
    oldBED = args.bfile + '.BED'
    makesure(bim2bed(args.bfile+'.bim', oldBED),
             'map->bed succ')

    # liftover to new BED file
    newBED = args.prefix + '.BED'
    unlifted = args.prefix + '.unlifted'
    makesure(liftBed(oldBED, newBED, unlifted, args.chainFile, args.liftOverPath),
             'liftOver succ')

    # extract lifted variants
    temp = tempfile.NamedTemporaryFile('w+t', delete=False)
    tmp_bfile = temp.name 
    makesure(extractBfile(args.bfile, newBED, tmp_bfile),
             'Extract Bfile succ')
    
    # update bim file with new position
    makesure(updatebim(tmp_bfile+'.bim', newBed),
             'update bim succ')

    # dry run PLINK to order bim file
    basic_fun.subprocess_popen(['plink', '--bfile', tmp_bfile, '--memory', args.memory, 
        '--make-bed', '--out', args.prefix])

    if args.datFile:
        newDat = args.prefix + '.dat'
        makesure(liftDat(args.datFile, newDat),
                 'liftDat succ')


    print("cleaning up BED files...")
    os.remove(newBed)
    os.remove(oldBed)

