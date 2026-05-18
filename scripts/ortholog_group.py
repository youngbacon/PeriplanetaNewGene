"""
    The script can generate ortholog/paralog mapping providing cactus alignments
    Be aware that halLiftover is not capabale of finding duplications of species-specific genes 
"""

import argparse
import sys

def bool_same_gene(coor_genei,coor_genej,cutoff):
    try:
        si,coori = coor_genei.split("|")
        sj,coorj = coor_genej.split("|")
    except ValueError:
        print(coor_genei,coor_genej)
        sys.exit()

    if coor_genei==coor_genej:
        return coor_genei

    else:
        flag = 0
        if si==sj and ":" in coor_genei and ":" in coor_genej:
            chromi,starti,endi = coori.split(":")
            strandi = chromi[-1]
            chromi  = chromi[0:-1]
 
            chromj,startj,endj = coorj.split(":")
            strandj = chromj[-1]
            chromj  = chromj[0:-1]
 
            if chromi==chromj and strandi==strandj:
 
                starti = int(starti)
                startj = int(startj)
                endi   = int(endi)  
                endj   = int(endj)  
      
                bound1 = max(starti,startj)-cutoff
                bound2 = min(endi,endj)    +cutoff
 
                if bound1 < bound2:
                    return "%s|%s%s:%d:%d"%(si,chromi,strandi,min(starti,startj),max(endi,endj))
    return False

class OrthoMap():
    def __init__(self,gene,cutoff=300):
        self.gene      = gene
        self.paralogs  = []
        self.orthologs = []
        self.has_dup   = False ## if this gene has a duplication?
        self.unique    = True  ## if this gene is species specific?
        self.cutoff    = cutoff

    def add_paralog(self,g):
        self.paralogs.append(g)

    def add_ortholog(self,g):
        self.orthologs.append(g)


class OrthologsMap():
    def __init__(self,fList,Genelst,cutoff=300):
        """ 
            Attributs:
             - self.key                  focal species
             - self.specific_genes      genes in focal species
             - self.orthomap            ortholog & paralog info of each gene

            note1: In the current state, I only compare dmel to other species

            note2: Ultimately, I want self.orthomap to contain all genes 
                   in all species
        """
        self.orthomap   = dict()
        self.cutoff     = cutoff

        species = []
        for i in open(fList,'r'):
            species.append(i.strip())
        self.key = species[0]
        self.species = species[1:]
        self.specific_genes = dict()
        lines = open(Genelst,"r")
        for line in lines:
            g = "%s|"%self.key + line.strip()
            self.specific_genes[g] = 1
            self.orthomap[g] = OrthoMap(g)
        lines.close()

        #self.uniqueGene = dict()

    def pairwiseOrthologs(self,s1,s2,pair_dir):
        """ 1. Identify duplications in s1 relative to s2
            2. Identify one-to-one ortholog between s1 and s2
        """

        genesets = dict()
        ## first assign all ortholog information
        fmap  = "%s/cactus_orthologs_%s_vs_%s.txt"%(pair_dir,s1,s2)
        lines = open(fmap,"r")
        for line in lines:
            linesplit = line[:-1].split()
            g1 = linesplit[0]
            genes_matched = linesplit[1:]

            if g1 not in genesets:
                genesets[g1] = 1 

            if g1 not in self.orthomap:
                self.orthomap[g1] = OrthoMap(g1)
                self.orthomap[g1].unique  = False

            for g2 in genes_matched:
                if g2 not in genesets:
                    genesets[g2] = 1

                ## add g2 to g1's orthologs
                if g2 not in self.orthomap[g1].orthologs:
                    self.orthomap[g1].add_ortholog(g2)

                ## also add g1 to g2's orthologs
                if g2 not in self.orthomap:
                    self.orthomap[g2] = OrthoMap(g2)
                    self.orthomap[g2].unique  = False
                    self.orthomap[g2].add_ortholog(g1)
                        
                else:
                    if g1 not in self.orthomap[g2].orthologs:
                        self.orthomap[g2].add_ortholog(g1)

            ## process duplications
            n_matched = len(genes_matched)
            if n_matched > 1:
                for i in range(n_matched-1):
                    for j in range(i+1,n_matched):
                        g1 = genes_matched[i]
                        g2 = genes_matched[j]

                        self.orthomap[g1].has_dup = True
                        self.orthomap[g2].has_dup = True

                        ## add g2 to g1's paralogs
                        if g2 not in self.orthomap[g1].paralogs:
                            self.orthomap[g1].add_paralog(g2)

                            ## in the meantime, add each of g1's paraplogs to
                            ## g2, and add g2 to each of g1's paralogs's paralog
                            for g3 in self.orthomap[g1].paralogs:
                                if g3 not in self.orthomap[g2].paralogs:
                                    self.orthomap[g2].add_paralog(g3)
                                if g2 not in self.orthomap[g3].paralogs:
                                    self.orthomap[g3].add_paralog(g2)


                        ## add g1 to g2's paralogs
                        if g1 not in self.orthomap[g2].paralogs:
                            self.orthomap[g2].add_paralog(g1)

                            ## in the meantime, add each of g2's paraplogs to
                            ## g1, and add g1 to each of g2's paralogs's paralog
                            for g3 in self.orthomap[g2].paralogs:
                                if g3 not in self.orthomap[g1].paralogs:
                                    self.orthomap[g1].add_paralog(g3)
                                if g1 not in self.orthomap[g3].paralogs:
                                    self.orthomap[g3].add_paralog(g1)
        lines.close()

    def update(self,pair_dir):
        """ update ortholog/paralog relations, all collapse onto focal species
        """
        ns  = len(self.species)
        for i in range(0,ns):
            s1 = self.key
            s2 = self.species[i]
            self.pairwiseOrthologs(s1,s2,pair_dir)

        groups    = dict()
        for g in self.specific_genes:
            groups[g] = self.orthomap[g].paralogs
            for g1 in self.orthomap[g].paralogs:
                if g1 in self.orthomap:
                    for g2 in self.orthomap[g1].paralogs:
                        if g2 not in self.orthomap[g].paralogs:
                            self.orthomap[g].paralogs.append(g2)
                        if g2 not in groups[g]:
                            groups[g].append(g2)
                    for g2 in self.orthomap[g1].orthologs:
                        if g2 not in self.orthomap[g].orthologs:
                            self.orthomap[g].orthologs.append(g2)

            for g1 in self.orthomap[g].orthologs:
                flag = 0
                for g0 in groups[g]:
                    combined = bool_same_gene(g0,g1,self.cutoff)
                    if combined:
                        groups[g].remove(g0)
                        groups[g].append(combined)
                        flag = 1
                if flag==0:
                    groups[g].append(g1)
                         
                for g2 in self.orthomap[g1].paralogs:
                    if g2 not in groups[g]:
                        flag = 0
                        for g0 in groups[g]:
                            combined = bool_same_gene(g0,g2,self.cutoff)
                            if combined:
                                groups[g].remove(g0)
                                groups[g].append(combined)
                                flag = 1
                        if flag==0:
                            groups[g].append(g2)

        ncopy = dict()
        groups_rm = []     ## duplications are included in only one group
        genes_rm  = dict() ## add a duplicate gene to genes_rm, and it will not
                           ## be considered again, this is to get groups_rm
        for s in self.species:
            ncopy[s] = 0
        idx = 0
        for g in self.specific_genes:
            if g not in genes_rm:
                if not self.orthomap[g].has_dup:
                    ncopy[self.key] = 1
                for g1 in groups[g]:
                    s = g1[0:4]
                    ncopy[s] += 1

                var = [ncopy[s] for s in self.species]
                if self.orthomap[g].has_dup:
                    groups_rm += [groups[g]]
                else:
                    groups_rm += [[g]+groups[g]]
                idx += 1
                    
                for g1 in self.orthomap[g].paralogs:
                    genes_rm[g1] = 1

                for s in self.species:
                    ncopy[s] = 0
        
        fout = open('cactus_ortholog_group.txt', 'w')
        for grp in groups_rm:
            newgrp = [s for s in grp]
            n_grp  = len(grp)
            for i in range(n_grp-1):
                for j in range(i+1,n_grp):
                    g1 = grp[i]
                    g2 = grp[j]
                    combined = bool_same_gene(g1,g2,self.cutoff)
                    if combined:
                        newgrp[i] = combined
                        newgrp[j] = ""
            newgrp = [s for s in newgrp if s!=""]
            
            for g in newgrp:
                s = g[0:4]
                ncopy[s] += 1
            
            var = [ncopy[s] for s in self.species]
            fout.write(" ".join(newgrp)+'\n')
                    
            for s in self.species:
                ncopy[s] = 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="bedAggregated")
    parser.add_argument("-f", "--fList",
                        help="Input file name")
    parser.add_argument("-p", "--Path",
                        help="Input tree file name")
    parser.add_argument("-g", "--Genelst",
                        help="Input key species name")
    # -------------------------------------------------------------------->>>>>
    # load the paramters
    # -------------------------------------------------------------------->>>>>
    ARGS = parser.parse_args()

    fList=ARGS.fList
    Path = ARGS.Path
    Genelst = ARGS.Genelst
    
    ortho = OrthologsMap(fList,Genelst)
    ortho.update(Path) 
