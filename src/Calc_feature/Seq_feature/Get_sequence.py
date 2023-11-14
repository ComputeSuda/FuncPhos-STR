import numpy as np
import pandas as pd 
import os
import sys


standard_amino_acids = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 
						'GLY': 'G', 'GLN': 'Q', 'GLU': 'E', 'HIS': 'H', 'ILE': 'I',
						'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PRO': 'P', 'PHE': 'F',
						'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'TRP': 'W', 'VAL': 'V'}


def get_pdb_site(pdb_file):
	site_no = []
	sequence = ''

	with open(pdb_file, 'r') as f:		
		for line in f:
			if line[:4] != 'ATOM':
				continue
			line_list = line.strip().split()
			sitenum = line[22:28].replace(' ', '')
			residue = line[17:20].replace(' ', '')  # THR PRO ...

			if residue not in standard_amino_acids.keys():
				print('错误：pdb文件找到了非标准氨基酸')
				print(sitenum)
				sys.exit()

			if sitenum not in site_no:
				site_no.append(sitenum)
				sequence += standard_amino_acids[residue]

		

	return sequence



if __name__ == '__main__':
	error_pdb=[]
	pdb_path=r'02structure\03del_low_pldd_acc'
	out_pdb_path = r'02structure\04del_low_pldd_acc_seq'
	pdb_list=os.listdir(pdb_path)
	for i in pdb_list:
		# print(i)
		pdb_file=os.path.join(pdb_path,i)
		# chain=i[5]
		pdb_fasta_name = i[:-4] + '.txt'
		fasta = os.path.join(out_pdb_path, pdb_fasta_name)
		try:
			sequence = get_pdb_site(pdb_file)
			with open(fasta,'w') as f:
				f.write('>'+i[:-4]+'_del_pldd'+'\n')
				f.write(sequence)
		except:
			error_pdb.append(i)
	print(error_pdb)









