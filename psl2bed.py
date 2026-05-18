# _*_ coding: UTF-8 _*_

# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    Author: YANG PEIWEN

    Version-01:
        2021-8 map seq
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
                        help="Input file name")
    parser.add_argument("-o", "--Output",
                        help="Output file name")
    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()
    full_cmd = """python pre_snpgenie.py\\
    -i {Input} \\
    -o {Output}\\
    """.format(
        Input=ARGS.Input,
        Output=ARGS.Output,
    )

    Input = ARGS.Input
    Output = ARGS.Output
    
    f2 = open(Output,"w")
    
    
    for line in open(Input, "r"):
        line = line[:-1].split('\t')
        if line[9] == '+-':
            if '+' in line[0]:
                line[0] = line[0][:-1]+'-'
            else:
                line[0] = line[0][:-1]+'+'
            f2.write(line[14]+'\t'+line[16]+'\t'+line[17]+'\t'+line[0]+'\n')
        else:
            f2.write(line[14]+'\t'+line[16]+'\t'+line[17]+'\t'+line[0]+'\n')
    
        
    
    
    f2.close()
    
    
# 2021-6-28















