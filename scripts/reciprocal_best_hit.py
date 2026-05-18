# _*_ coding: UTF-8 _*_

# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    Author: YANG PEIWEN

    Version-01:
        2021-6 map orthologs
    """
# Version information END ----------------------------------------------------

import argparse
import math
import numpy as np
import scipy.stats as stats
###############################################################################
# read parameters
###############################################################################


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="run snpgenie")
    parser.add_argument("-i", "--Input",
                        help="Input file name")
    parser.add_argument("-p", "--Path",
                        help="Input file name")
    parser.add_argument("-o", "--Output",
                        help="Output file name")
    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()

    Input=ARGS.Input
    Path=ARGS.Path
    Output=ARGS.Output
    
    f2 = open(Output, "w")
    
    species = []
    for line in open(Input, "r"):
        species.append(line.strip())
        
    key = species[0]
    d = {}
    dmax = {}
    for i in species[1:]:
        d = {}
        dmax = {}
        for line in open(Path+key+'_'+i+'.blastp', "r"):
            if float(line[:-1].split()[-2]) > 1e-5:
                continue
            g1,g2 = line[:-1].split()[0:2]
            score = line[:-1].split()[-1]
            if g1 not in d.keys():
                dmax[g1] = float(score)
                d[g1] = g2
            elif float(score) > dmax[g1]:
                dmax[g1] = float(score)
                d[g1] = g2
                
            if g2 not in d.keys():
                dmax[g2] = float(score)
                d[g2] = g1
            elif float(score) > dmax[g2]:
                dmax[g2] = float(score)
                d[g2] = g1
    
        for i in d:
            for j in d:
                if key in i and key not in j:
                    if d[i] == j and d[j] == i:
                        f2.write(i+'\t'+j+'\n')
            
    
    
    f2.close()

# 2020-1-3










