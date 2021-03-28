#!/home/local/users/jw/anaconda2/bin/python
#Author: jywang	explorerwjy@gmail.com

#========================================================================================================
# select_gnomad_var.py
#========================================================================================================

import argparse
import gzip as gz
import pandas as pd
import csv
import pysam

LGD_Category = set(["splice_acceptor_variant", "splice_donor_variant", "stop_gained", 
				"stop_lost", "start_lost", "frameshift_variant"])

class select_gnomad_var:
	def __init__(self, args):
		self.gnomad = args.gnomad
		self.out_SYN = "gnomad.mut.1e-3.syn.bed"
		self.out_LGD = "gnomad.mut.1e-3.lgd.bed"
		self.out_MIS = "gnomad.mut.1e-3.mis.bed"
		self.AF_cut = 1e-3
	def select_var(self):
		hand = gz.open(self.gnomad, 'rt')
		outfil_syn = csv.writer(open(self.out_SYN, 'wt'), delimiter="\t")
		outfil_lgd = csv.writer(open(self.out_LGD, 'wt'), delimiter="\t")
		outfil_mis = csv.writer(open(self.out_MIS, 'wt'), delimiter="\t")
		head = ["#Chr", "Start", "End", "Gene", "AC", "Cons"]
		outfil_syn.writerow(head)
		outfil_lgd.writerow(head)
		outfil_mis.writerow(head)
		for l in hand:
			if l.startswith("##INFO=<ID=vep"):
				CSQ_header = l.strip().split("Format: ")[1].rstrip('>\"').split("|")
			elif l.startswith("#"):
				continue
			else:
				try:
					llist = l.strip().split("\t")
					Chr, Pos, Ref, Alts = llist[0], llist[1], llist[3], llist[4]
					if llist[6] != "PASS":
						continue
					INFO = VCF_INFO_Parser(llist[7])
					Allele_CSQ_dict = VCF_VEP_CSQ_Parser(Ref, Alts, CSQ_header, INFO["vep"])
					for i, alt in enumerate(Alts.split(',')):
						AC = int(INFO["AC"].split(",")[i])
						AF = float(INFO["AF"].split(",")[i])
						#print(AF, self.AF_cut)
						if AF > self.AF_cut:
							continue
						vep = Allele_CSQ_dict[alt][0]
						cons = vep["Consequence"]
						Gene_Symbol = vep["SYMBOL"]
						Start, End = int(Pos), int(Pos) + len(alt)
						Var = "{}-{}-{}".format(Pos, Ref, alt)
						if len(set(cons).intersection(LGD_Category))>= 1:
							outfil_lgd.writerow([Chr, Start, End, Var, Gene_Symbol, AC, ",".join(cons)]) #LGD
						if "synonymous_variant" in set(cons):
							outfil_syn.writerow([Chr, Start, End, Var, Gene_Symbol, AC, ",".join(cons)]) #SYN
						if "missense_variant" in set(cons):
							outfil_mis.writerow([Chr, Start, End, Var, Gene_Symbol, AC, ",".join(cons)]) #MIS
				except:
					print(l)
					

def VCF_INFO_Parser(info_string):
	infolist = info_string.split(';')
	infodict = {}
	for kv in infolist:
		kv = kv.split('=')
		if len(kv) == 2:
			k,v = kv
			infodict[k] = v
	return infodict

def VCF_VEP_CSQ_Parser(Ref, Alts, csq_head, csq_string):
	Alts = Alts.split(",")
	if len(list(set([x[0] for x in Alts])))==1 and Ref[0] == list(set([x[0] for x in Alts]))[0]:
		_Ref = Ref[1:] if len(Ref[1:]) >0 else "-"
		_Alts = [Alt[1:] if len(Alt[1:]) >0 else "-" for Alt in Alts]
	else:
		_Alts = Alts
	res = {}
	csqs = csq_string.split(",")
	csqs = [dict(zip(csq_head, vep.split("|"))) for vep in csqs]
	for i, Alt in enumerate(Alts):
		res[Alt] = []
		for j, csq in enumerate(csqs):
			if csq["Allele"] == _Alts[i]:
				csq["Consequence"] = csq["Consequence"] .split("&")
				res[Alt].append(csq)
	return res

def gtf_info_parser(info):
    res = {}
    for term in info.split(";"):
        if term == "":
            continue
        #print(">",term)
        key,v = term.split()
        v = v.strip('"')
        res[key]=v
    return res

def GetOptions():
	parser = argparse.ArgumentParser()
	parser.add_argument('--gnomad', default="/home/local/users/jw/resources/gnomAD/exome_2.1/gnomad.exomes.r2.1.1.sites.vcf.gz", type=str, help = 'gnmoAD file')
	args = parser.parse_args()
	
	return args

def main():
	args = GetOptions()
	ins = select_gnomad_var(args)
	ins.select_var()	

	return

if __name__=='__main__':
	main()
