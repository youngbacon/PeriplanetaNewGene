# _*_ coding: UTF-8 _*_

# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    Author: YANG PEIWEN

    Version-01:
        2024-7 map seq
    """
# Version information END ----------------------------------------------------

import argparse
import math
import scipy.stats as stats

###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="map seq")
    parser.add_argument("-i", "--Input",
                        help="Input  file name")
    parser.add_argument("-o", "--Output",
                        help="Output file name")
    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()

    Input = ARGS.Input
    Output = ARGS.Output
    
    f2 = open(Output,"w")
    
    for line  in open(Input, "r"):
        line = line.strip().split()
        if line[1] in  line[0]:
            l = line[0].split('_')
            f2.write('\t'.join(l)+'\n')
    
    
    f2.close()
    
    
# 2021-6-28























