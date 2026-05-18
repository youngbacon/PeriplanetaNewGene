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

class Stack:
    def __init__(self):
        self.items = []
    def isEmpty(self): 
        return self.items == []
    def push(self, item):
        self.items.append(item)
    def pop(self):
        return self.items.pop()
    def peek(self):
        return self.items[len(self.items)-1]
    def size(self):
        return len(self.items)

class Tree():
    def __init__(self):
        self.items = [[],[]]
        self.leaf = False
        self.parent = None
    def getChild(self, n):
        return self.items[n]
    def getParent(self):
        return self.parent
    def insert(self, n = 0, N = None):
        if N == None:
            N = Tree()
        self.items[n] = N
        N.setParent(self)
    def setParent(self,N):
        self.parent = N
    def setVal(self, Val):
        self.leaf = True
        self.items = Val
    def getVal(self):
        if self.leaf == True:
            return self.items
        else:
            return self.leaf
    def isEmpty(self):
        if self.items == [[],[]]:
            return True
        else:
            return False
    def getPos(self,N):
        if N == self.items[0]:
            return 0
        elif N == self.items[1]:
            return 1
        else:
            return False
        
def BuildRooted(l):
    e = Tree()
    p = Stack()
    currentTree = e
    for i in l:
        if i == '(':
            p.push(currentTree)
            currentTree.insert(0)
            currentTree.insert(1)
            currentTree = currentTree.getChild(0)
            node = ''
        elif i == ')':
            if node != '':
                currentTree.setVal(node)
                node = ''
            parent = p.pop()
            currentTree = parent
        elif i == ',':
            if node != '':
                currentTree.setVal(node)
                node = ''
            parent = p.peek()
            if parent.getChild(1).isEmpty():
                currentTree = parent.getChild(1)
            else:
                parent = p.pop()
                newroot = Tree()
                newroot.insert(0,parent)
                newroot.insert(1)
                p.push(newroot)
                currentTree = newroot.getChild(1)
        elif i == ';':
            break
        else:
            node += i
    return currentTree

def fetch(tree, v):
    if tree.getVal() == v:
        return tree
    elif not tree.getVal():
        return fetch(tree.getChild(0), v) or fetch(tree.getChild(1), v)
    else:
        return False


def parseTree(T, Key):
    tree = BuildRooted(T)
    TreeD = {}
    species = T[:-1].split(',')
    
    kn = fetch(tree, Key)
    for i in species:
        i = i.replace('(','').replace(')','').replace(';','')
        currentNode = kn
        c = 1
        while currentNode:
            if fetch(currentNode, i):
                break
            else:
                currentNode = currentNode.getParent()
                c += 1
        TreeD[i] = c
    
    return TreeD

def ReadRBH(RBH):
    rbh = {}
    for line in open(RBH):
        line = line.strip().split()
        if line[0] not in rbh:
            rbh[line[0]] = [line[1]]
        else:
            rbh[line[0]].append(line[1])
    return rbh

def ReadOFG(OFG):
    c = {}
    for line in open(OFG, "r"):
        line = line[:-1].split('\t')
        tag = line[1][-1]
        for i in line[2:]:
            if 'Pame' in i:
                gene = i.split(',')[0]
                c[gene] = tag
    
    return c

def ReadBlastp(fblastp):
    """ parse similar sequence to potential paralog pairs
    """
    para_pairs = dict()
    lines = open(fblastp,"r")
    for line in lines:
        g1,g2 = line.split()[0:2]
        if g1 == g2:
            continue
        if g1 not in para_pairs:
            para_pairs[g1] = [g2]
        else:
            if g2 not in para_pairs[g1]:
                para_pairs[g1].append(g2)

        if g2 not in para_pairs:
            para_pairs[g2] = [g1]
        else:
            if g1 not in para_pairs[g2]:
                para_pairs[g2].append(g1)
    lines.close()
    
    return para_pairs


###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="map seq")
    parser.add_argument("-i", "--Input",
                        help="Input cactus group file name")
    parser.add_argument("-t", "--Tfile",
                        help="Tree file name")
    parser.add_argument("-s", "--Seq",
                        help="Orthoseq file name")
    parser.add_argument("-r", "--Reciprocal",
                        help="Reciprocal best hit file name")
    parser.add_argument("-b", "--blastp",
                        help="Self-blastp file name")
    parser.add_argument("-o", "--Output",
                        help="Output file name")
    parser.add_argument("-k", "--KeyS",
                        help="Key Species name, seperated by comma; for example: Pame,Pjap")
    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()

    Input = ARGS.Input
    Seq = ARGS.Seq
    Tfile = ARGS.Tfile
    Reciprocal = ARGS.Reciprocal
    blastp = ARGS.blastp
    Output = ARGS.Output
    KeyS = ARGS.KeyS
    
    f2 = open(Output,"w")
    
    T = open(Tfile, 'r').readlines()[0]
    s1 = KeyS.split(',')[0]
    s2 = KeyS.split(',')[1]
    TreeDis = parseTree(T, s1)
    TreeCut = TreeDis[s2]
    
    para_pairs = ReadBlastp(blastp)
    rbh = ReadRBH(Reciprocal)
    
    paras = []
    agf = []
    d = {}
    wlst = []
    #Check gene duplication or annotation error/gene fission/fusion
    for line in open(Input, "r"):
        line = line[:-1].split()
        genes = [i.split('|')[1] for i in line]
        
        keygenes = []
        for i in genes:
            if s1 in i:
                keygenes.append(i)
        
        if len(keygenes) > 1:
            for i in keygenes:
                for j in keygenes:
                    if i != j and i in para_pairs:
                        if j in para_pairs[i]:
                            paras.append(i)
                            paras.append(j)
            for i in keygenes:
                if i not in paras:
                    agf.append(i)
        
        for i in keygenes:
            d[i] = []
            for j in genes:
                if ':' not in j:
                    d[i].append(j[0:4])
        wlst += genes
    
    paras = set(paras)
    
    for i in rbh:
        for j in rbh[i]:
            if j not in wlst:
                d[i].append(j[0:4])
                    
    for line in open(Seq, "r"):
        line = line[:-1].split()
        d[line[0]].append(line[1][0:4])
    
    # catagory genes
    for i in d:
        clade = max([TreeDis[j] for j in d[i]])
        if i in agf:
            f2.write(i+'\tAnnotationErrorFusionFission\t'+str(clade)+'\n')
        elif i in paras:
            f2.write(i+'\tDup_Candidate\t'+str(clade)+'\n')
        elif clade > TreeCut:
            f2.write(i+'\tConserved\t'+str(clade)+'\n')
        elif i not in para_pairs:
            f2.write(i+'\tOrphan\t'+str(clade)+'\n')
        else:
            f2.write(i+'\tDup_Candi_2\t'+str(clade)+'\n')
    
    
    
    
    f2.close()
    
    
# 2021-6-28






















