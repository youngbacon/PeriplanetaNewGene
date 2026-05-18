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
                        help="Input catagory file name")
    parser.add_argument("-o", "--Output",
                        help="Output file name")
    parser.add_argument("-a", "--Anno", default='pep',
                        help="Output file name")
    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()

    Input = ARGS.Input
    Output = ARGS.Output
    Anno=ARGS.Anno
    
    f2 = open(Output,"w")
    
    for line in open(Input, 'r'):
        line = line[:-1].split()
        check = []
        for i in line:
            if 'Pame' in i:
                check.append(i.split('|')[1])
            else:
                if ':' in i:
                    for j in check:
                        g = i.split('|')[0]
                        scff, stt, end = i.split('|')[1].split(':')
                        strand = scff[-1]
                        scff = scff[:-1]
                        f2.write("/data/yangpeiwen/application/seqkit grep -p " +scff+" /data/yangpeiwen/SP19/WGS/"+g+".chromosome.fa |")
                        f2.write("/data/yangpeiwen/application/seqkit subseq -r "+stt+":"+end+' > tmp.genome \n')
                        f2.write("/data/yangpeiwen/application/seqkit grep -p " +j+ " /data/yangpeiwen/SP19/03_New_gene/01_WGA/Pame.protein.fa > tmp.pep \n")
                        if(strand == '+'):
                            f2.write("/data/yangpeiwen/application/wise2.4.1/src/bin/genewise tmp.pep tmp.genome -pseudo -genes -"+Anno+" -quiet -tfor -gff > check_seq/"+j+'_'+g+" \n")
                        else:
                            f2.write("/data/yangpeiwen/application/wise2.4.1/src/bin/genewise tmp.pep tmp.genome -pseudo -genes -"+Anno+" -quiet -trev -gff > check_seq/"+j+'_'+g+" \n")
                        
    
    f2.close()
    
    
# 2021-6-28





















