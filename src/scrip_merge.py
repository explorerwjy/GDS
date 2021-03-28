from GDS import *
import pandas as pd
Gencode_exon = pd.read_csv("/Users/jiayao/Work/BrainDisorders/data/GTEx/gencode.v19.genes.v7.patched_contigs.exons.txt",
                          delimiter="\t", index_col="exon_id")
GeneExon = pd.read_csv("/Users/jiayao/Work/BrainDisorders/data/GTEx/GENE_Exon.map", index_col="Name")
for i, row in Gencode_exon.iterrows():
    Gencode_exon.loc[i, "gene"] = GeneExon.loc[i, "Description"]
Gencode_exon.to_csv("test.tsv", delimiter="\t")

