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
    parser.add_argument("-e", "--expressed",
                        help="Input expression file name")
    parser.add_argument("-u", "--uniprot",
                        help="Input uniprot blastp file name")
    parser.add_argument("-o", "--Output",
                        help="Output file name")
    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()

    Input = ARGS.Input
    E = ARGS.expressed
    U = ARGS.uniprot
    Output = ARGS.Output
    
    f2 = open(Output,"w")
    tissue = ['AG','AT','Br','Fcg','M10','M1','M4','M7','Mcg','OV','Testis']
    eg = {}
    d = {}
    for i in tissue:
        d[i] = []
    for line in open(E):
        if line[0] != 'P':
            line = line.strip().split()
            for i in range(1,len(line)):
                for j in tissue:
                    if j in line[i]:
                        d[j].append(i)
        else:
            line = line.strip().split()
            flag = 0
            for i in d:
                c = 0
                for j in d[i]:
                    c += float(line[j])
                if c/len(d[i]) > flag:
                    flag = c/len(d[i]) 
            
            if flag > 1:
                eg[line[0]] = flag
    
    morehit = []
    for line in open(U, "r"):
        line = line[:-1].split()
        morehit.append(line[0])
    morehit = set(morehit)
    
    for line in open(Input, "r"):
        line = line.split()
        if line[0] in morehit:
            line[1] = 'Dup_Candi_3'
        line = '\t'.join(line)
        if line[0:9] in eg:
            f2.write(line+'\t'+str(eg[line[0:9]])+'\n')
    
    
    f2.close()
    
    
# 2021-6-28






















