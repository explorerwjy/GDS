import numpy as np
import pysam 
import matplotlib.pyplot as plt
import pandas as pd
import csv
import statsmodels.api as sm

########################################################################################################
# Mutation Distribution Across gene
########################################################################################################
class Variant():
    def __init__(self, ID, Chr, Pos, AC):
        self.ID = str(Chr) + "-" + ID
        self.Chr = Chr
        self.Pos = int(Pos)
        self.AC = int(AC)

class ExonVar():
    def __init__(self, ExonID, Start, End):
        self.ID = ExonID
        self.Start = int(Start) - 3
        self.End = int(End) + 3
        self.Vars = []
    def addVar(self, Var):
        self.Vars.append(Var)
    def countVar(self):
        return sum([var.AC for var in self.Vars])

class GeneExonVar:
    def __init__(self, gene, df_exon):
        self.Gene = gene
        self.Strand = df_exon["strand"].values[0]
        self.Exons = []
        for i, row in df_exon.iterrows():
            Exon = ExonVar(i, row["start_pos"], row["end_pos"])
            self.Exons.append(Exon)
        if self.Strand == "-":
            self.Exons = self.Exons[::-1]
    def ShowExons(self, withVar=False):
        for exon in self.Exons:
            print(exon.Start, exon.End, exon.countVar())
    def addVars(self, df_Var): #Exon and Variants must be sorted
        idx_exon = 0
        for i, row in df_Var.iterrows():
            Var = Variant(row["VAR"], row["#Chr"], row["Start"], row["AC"])
            SKIP_VAR = False
            while (idx_exon < len(self.Exons)) and (not SKIP_VAR):
                # ---||||||------||||||--
                #  x
                if Var.Pos < self.Exons[idx_exon].Start:
                    idx_exon += 1
                # ---||||||------||||||--
                #            x
                elif Var.Pos > self.Exons[-1].End: # Intergenic
                    break
                elif Var.Pos > self.Exons[idx_exon].End and Var.Pos < self.Exons[idx_exon+1].Start:
                    #print("Intron")
                    SKIP_VAR = True
                    break
                elif Var.Pos >= self.Exons[idx_exon].Start and Var.Pos <= self.Exons[idx_exon].End:
                    self.Exons[idx_exon].addVar(Var)
                    break
                else:
                    idx_exon += 1
            #if idx_exon >= len(self.Exons):
            #    self.ShowExons()
            #    print(Var.ID)
    def Percentile(self, percent_step=10):
        Cum_Len = 0
        POS, ACs = [], []
        for exon in self.Exons:
            for var in exon.Vars:
                rel_pos = Cum_Len + var.Pos - exon.Start
                POS.append(rel_pos)
                ACs.append(var.AC)
            Cum_Len += exon.End - exon.Start
        POS = np.array(POS)
        ACs = np.array(ACs)
        Percents = POS/Cum_Len
        Mut_Den = ACs/sum(ACs)
        #print(Cum_Len, Percents, Mut_Den)
        return Percents, Mut_Den

def Percentile(Percents, Mut_Den, l):
    #l = np.arange(0, 1, percent_step) + percent_step
    percent_idx = 0
    last_v = 0
    res = np.zeros(len(l))
    for i, v in enumerate(l[:-1]):
        if percent_idx >= len(Percents):
            break
        #print(i, percent_idx)
        while percent_idx < len(Percents):
            #print(Percents[percent_idx], [last_v, v])
            if (Percents[percent_idx] >= last_v) and (Percents[percent_idx] <= v):
                res[i] += Mut_Den[percent_idx]
                percent_idx += 1
                # stay in loop, see if next percent still in v_range
            elif Percents[percent_idx] < v:
                #continue
                break
            elif Percents[percent_idx] > v:
                #percent_idx += 1
                break
            elif Percents[percent_idx] < last_v:
                #print("Error I")
                break
        last_v = v
    return res
            

########################################################################################################
# GDS 
########################################################################################################
def poisson_reg(X_train, y_train):
    res = sm.GLM(y_train, X_train, family=sm.families.Poisson()).fit()
    return res

def GDS_one_gene_tissue(gene, Tissue, Dat, rel_exp):
    exons = {}
    for exon in Dat[gene].Exons:
        exons[exon.ID] = exon
    tmp_x = rel_exp[rel_exp["Description"]==gene]
    Tissues = tmp_x.columns.values[1:-1]
    for i, row in tmp_x.iterrows():
        exon_id = row["Name"]
        tmp_x.loc[i, "intercept"] = 1
        tmp_x.loc[i, "N_mut"] = np.sum([var.AC for var in exons[exon_id].Vars])
    x_train = tmp_x[Tissue].values
    y_train = tmp_x["N_mut"].values
    print(x_train)
    plt.scatter(x_train, y_train)


def GDS_one_gene(Dat, gene, rel_exp):
    exons = {}
    for exon in Dat[gene].Exons:
        exons[exon.ID] = exon
    tmp_x = rel_exp[rel_exp["Description"]==gene]
    Tissues = tmp_x.columns.values[1:-1]
    for i, row in tmp_x.iterrows():
        exon_id = row["Name"]
        tmp_x.loc[i, "intercept"] = 1
        tmp_x.loc[i, "N_mut"] = np.sum([var.AC for var in exons[exon_id].Vars])
    for Tissue in Tissues:
        #print(Tissue)
        x_train = tmp_x[["intercept", Tissue]]
        #x_train[Tissue] = np.exp(x_train[Tissue])
        #print(x_train[Tissue])
        y_train = tmp_x["N_mut"]
        res = poisson_reg(x_train, y_train)
        alpha, beta = res.params["intercept"], res.params[Tissue]
        p_alpha, p_beta = res.pvalues["intercept"], res.pvalues[Tissue]
        print(Tissue, beta, p_beta) 





