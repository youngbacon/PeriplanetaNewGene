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
    
    for file in open(Input, "r"):
        flag = False
        gene = file.strip().split('/')[1]
        for line in open(file.strip(), "r"):
            if line[0] == "/":
                flag = False
            elif line[0] == '>':
                flag = True
            
            if flag and line[0] == '>':
                f2.write('\n>'+gene+'\n')
            elif flag:
                f2.write(line)
    
    
    
    f2.close()
    
    
# 2021-6-28






















