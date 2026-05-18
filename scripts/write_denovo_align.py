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
    parser.add_argument("-l", "--Lst",
                        help="Input catagory file name")
    parser.add_argument("-o", "--Output",
                        help="Output file name")
    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()

    Input = ARGS.Input
    Lst = ARGS.Lst
    Output = ARGS.Output
    
    f2 = open(Output,"w")
    
    gene = {}
    for line in open(Input, 'r'):
        gene[line[:-1]] = {}
    
    for i in open(Lst, 'r'):
        s = i[:-1]
        for line in open('../01_WGA/bed_halLifted/Pame_'+s+'_cds.bed', "r"):
            line = line.strip().split()
            if line[3] in gene:
                if s not in gene[line[3]]:
                    gene[line[3]][s] = {}
                    gene[line[3]][s]['line'] = []
                gene[line[3]][s]['strand'] = line[-1]
                gene[line[3]][s]['line'].append(line)
                gene[line[3]][s]['scff'] = line[0]
    
    strand= {}
    for line in open('../01_WGA/bed_original/Pame_cds.bed', "r"):
        line = line.strip().split()
        if line[3] in gene:
            strand[line[3]]= line[-1]

    for i in gene:
        f2.write("/data/yangpeiwen/application/seqkit grep -p " +i+ " /data/yangpeiwen/SP19/WGS/Pame.cds.fa > tmp1.genome \n")
        for s in gene[i]:
            f2.write("echo \">"+s+"\" >> tmp1.genome \n")
            scff = gene[i][s]['scff']
            l1 = []
            l2 = []
            for line in gene[i][s]['line']:
                l1.append(int(line[1]))
                l2.append(int(line[2]))
            if gene[i][s]['strand'] == '+':
                l1.sort()
                l2.sort()
                for n in range(len(l1)):
                    f2.write("/data/yangpeiwen/application/seqkit grep -p " +scff+" /data/yangpeiwen/SP19/WGS/"+s+".chromosome.fa | ")
                    if strand[i] == '+':
                        f2.write("/data/yangpeiwen/application/seqkit subseq -r "+str(l1[n])+":"+str(l2[n]) )
                    else:
                        f2.write("/data/yangpeiwen/application/seqkit subseq -r "+str(l1[n]+1)+":"+str(l2[n]+1) )
                    f2.write(" | grep -v \'>\' >> tmp1.genome \n")
            elif gene[i][s]['strand'] == '-':
                l1.sort(reverse = True)
                l2.sort(reverse = True)
                for n in range(len(l1)):
                    f2.write("/data/yangpeiwen/application/seqkit grep -p " +scff+" /data/yangpeiwen/SP19/WGS/"+s+".chromosome.fa | ")
                    if  strand[i]  == '-':
                        f2.write("/data/yangpeiwen/application/seqkit subseq -r "+str(l1[n])+":"+str(l2[n]) )
                    else:
                        f2.write("/data/yangpeiwen/application/seqkit subseq -r "+str(l1[n]+1)+":"+str(l2[n]+1) )
                    f2.write(" | /data/yangpeiwen/application/seqkit seq -r -p | grep -v \'>\'>> tmp1.genome \n")
        f2.write('/data/yangpeiwen/application/MUSCLE/muscle -in tmp1.genome -out '+i+'.muscle \n')
        
    f2.close()
    
    
# 2021-6-28
























