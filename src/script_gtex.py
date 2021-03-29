#!/home/local/users/jw/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# script_gtex.py
# Calculate RPK Matrix
# calculate TPM Matrix
#========================================================================================================

import argparse
import csv
import gzip as gz
import numpy as np
import pandas as pd

class script_gtex:
	def __init__(self, args):
		pass
	def calculateRPK(self):
		f_exon_rc = "/home/local/users/jw/GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_exon_reads.csv.gz"
		reader = csv.reader(gz.open(f_exon_rc, 'rt'))
		writer = csv.writer(open("GTEx_v7_RNASeQCv1.1.8_exon_rpk.csv", 'wt'))
		df_gencode_exon = pd.read_csv("~/GTEx/gencode.v19.genes.v7.patched_contigs.exons.txt", 
				delimiter="\t", index_col="exon_id")
		head = next(reader)
		writer.writerow(head)
		for r in reader:
			#gene, ID = r[0], r[-1]
			exon = df_gencode_exon.loc[r[-1], :]
			exon_length = abs(int(exon["end_pos"]) - int(exon["start_pos"]) + 1)
			rpk_r = r[::]
			for i in range(1, len(r)-1, 1):
				rpk_r[i] = float(rpk_r[i]) * 1000 / exon_length
			writer.writerow(rpk_r)
	def calculateTPM(self):
		df_rpk = pd.read_csv("GTEx_v7_RNASeQCv1.1.8_exon_rpk.csv")
		writer = csv.writer(open("GTEx_v7_RNASeQCv1.1.8_exon_scale.csv", 'wt'))
		writer.writerow(["Sample", "Scale"])
		samples = df_rpk.columns.values[1:-1]
		for sample in samples:
			scale = df_rpk[sample].sum()
			scale = scale/1000000
			writer.writerow([sample, scale])
			df_rpk[sample] = df_rpk[sample] / scale
		df_rpk.to_csv("GTEx_v7_RNASeQCv1.1.8_exon_tpm.csv", index=False)
	def calculateTissueExpExon(self):
		reader = csv.reader(open("GTEx_v7_RNASeQCv1.1.8_exon_tpm.csv", 'rt'))
		writer = csv.writer(open("GTEx_v7_RNASeQCv1.1.8_exon_TissueExp.csv", 'wt'))
		GTEx_sample_attr = pd.read_csv("../dat/GTEx_v7_Annotations_SampleAttributesDS.txt", sep="\t")
		TissueList = sorted(list(set(GTEx_sample_attr["SMTS"].values)))
		TissueSampleDict = {}
		for Tissue in TissueList:
			SampleSet = set(GTEx_sample_attr[GTEx_sample_attr["SMTS"]==Tissue]["SAMPID"].values)
			TissueSampleDict[Tissue] = SampleSet
		head = next(reader)
		writer.writerow([head[0]] + TissueList + [head[-1]])
		for r in reader:
			dat = dict(zip(head, r))
			Tissue_Exp = []
			for Tissue in TissueList:
				SampleSet = TissueSampleDict[Tissue]
				explist = np.array([float(v) for (k,v) in dat.items() if k in SampleSet])
				explist = np.log(explist+1)
				logexp = np.mean(explist)
				Tissue_Exp.append(logexp)
			new_r = [r[0]] + Tissue_Exp + [r[-1]]
			writer.writerow(new_r)
	def calculateTissueExpGene(self):
		reader = csv.reader(open("/home/local/users/jw/GTEx/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct", 'rt'), delimiter="\t")
		writer = csv.writer(open("GTEx_v7_RNASeQCv1.1.8_gene_TissueExp.csv", 'wt'))
		GTEx_sample_attr = pd.read_csv("../dat/GTEx_v7_Annotations_SampleAttributesDS.txt", sep="\t")
		TissueList = sorted(list(set(GTEx_sample_attr["SMTS"].values)))
		TissueSampleDict = {}
		for Tissue in TissueList:
			SampleSet = set(GTEx_sample_attr[GTEx_sample_attr["SMTS"]==Tissue]["SAMPID"].values)
			TissueSampleDict[Tissue] = SampleSet
		head = next(reader)
		writer.writerow(head[0:2] + TissueList)
		for r in reader:
			dat = dict(zip(head, r))
			Tissue_Exp = []
			for Tissue in TissueList:
				SampleSet = TissueSampleDict[Tissue]
				explist = np.array([float(v) for (k,v) in dat.items() if k in SampleSet])
				if len(explist) == 0:
					print(Tissue, SampleSet)
				explist = np.log(explist+1)
				logexp = np.mean(explist)
				Tissue_Exp.append(logexp)
			new_r = r[0:2] + Tissue_Exp
			writer.writerow(new_r)
	def calculateRelativeExp(self):
		GTEx_Exon_TissueTPM = pd.read_csv("../dat/GTEx_v7_RNASeQCv1.1.8_exon_TissueExp.csv")
		GTEx_Gene_TissueTPM = pd.read_csv("../dat/GTEx_v7_RNASeQCv1.1.8_gene_TissueExp.csv", index_col="Name")
		writer = csv.writer(open("../dat/GTEx_v7_RNASeQCv1.1.8_exon_Tissue.RelExp.csv", "wt"))
		writer.writerow(GTEx_Exon_TissueTPM.columns.values)
		for i, row in GTEx_Exon_TissueTPM.iterrows():
			gene = row["Name"].split("_")[0]
			gene_exps = GTEx_Gene_TissueTPM.loc[gene, :].values[1:] 
			row_val = row.values[1:-1]
			rel_exp = row_val * 10 - gene_exps
			new_row = [row[0]] + list(rel_exp) + [row[-1]]
			writer.writerow(new_row)

def GetOptions():
	parser = argparse.ArgumentParser()
	#parser.add_argument('-f','--file', type=str, help = 'input file, 1. rc for rpk; 2. rpk for scaling and 3. rpk for tpm')
	#parser.add_argument('-s', '--scale', type=str, help = 'scaling factor file')
	parser.add_argument('-m','--method', type=str, help = 'method 1 for rpk; 2 for scaling factor and 3 for tpm')
	args = parser.parse_args()
	return args

def main():
	args = GetOptions()
	ins = script_gtex(args)
	if args.method == "1":
		print("Calculate RPK Matrix")
		ins.calculateRPK()
	elif args.method == "2":
		print("Calculate TPM Matrix")
		ins.calculateTPM()
	elif args.method == "3":
		ins.calculateTissueExpExon()
	elif args.method == "4":
		ins.calculateTissueExpGene()
	elif args.method == "5":
		ins.calculateRelativeExp()
	return

if __name__=='__main__':
	main()
