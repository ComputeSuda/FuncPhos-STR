import numpy as np
import pandas as pd 
import os
import sys
# from get_pdb_binding import get_pdb_binding_site

# 计算二级结构及ACC

def new_file(pdb_name,chain,name,data):
	dataframe=pd.DataFrame({name:data})
	csv_name=pdb_name+'_feature.csv'
	feature_path = r'dssp'
	csvname=os.path.join(feature_path,csv_name)
	dataframe.to_csv(csvname,index=False)

def add_content(pdb_name,chain,name, data):
	csv_name = pdb_name + '_feature.csv'
	feature_path =r'dssp'
	csvname = os.path.join(feature_path, csv_name)
	file = pd.read_csv(csvname)
	file[name] = data
	file.to_csv(csvname,index=False)


def run_dssp(pdb_file):
	pdb_name = pdb_file.split('/')[7][:-4]
	output_path = 'dssp_out/'
	if os.path.exists(output_path) == False:
		os.makedirs(output_path)
	output_file = output_path + pdb_name + '_dssp.out'
	os.system('dssp -i ' + pdb_file + ' -o ' + output_file)  # 在pdb文件上运行dssp



def get_second_structure_acc(pdb_file, chain):
	pdb_name = pdb_file[-10:-4]
	output_path = r'dssp_out/'
	dssp_out_file = output_path + pdb_name + '_dssp.out'

	sites = []
	with open(pdb_file, 'r') as f:
		for line in f:
			if line[:4] != 'ATOM':
				continue
			# if line[21] != chain:
			# 	continue
			line = line.strip()
			sitenum = line[22:28].replace(' ', '')

			if sitenum not in sites:
				sites.append(sitenum)
	
	site_ss = {}
	site_acc = {}
	# print(sites)
	ss=[str(0) for i in range(len(sites))]
	acc=[str(0) for i in range(len(sites))]
	# print(ss, len(ss))
	# print(acc, len(acc))
	with open(dssp_out_file, 'r') as f:
		start_flag = False
		pre_count = -1
		for line in f:
			line = line.rstrip()
			if line[2] != '#' and start_flag == False:
				continue
			if line[2] == '#' and start_flag == False:
				start_flag = True
				continue
			# if line[11] != chain:
			# 	continue
			site = line[5:11].replace(' ', '')  # 一级结构
			# print(site)
			s=line[16]
			ac=line[34:38].replace(' ','')
			if site in sites:
				count=sites.index(site)
				if s!=' ':
					ss[count]=s
				else:
					ss[count] = '0'
				acc[count]=ac
	# print(ss,len(ss))
	# print(acc,len(acc))



			#
			# if ss == ' ':
			# 	ss = '-'
			#
			# count = sites.index(site) + 1
			#
			# if count not in site_ss.keys():
			# 	site_ss[count] = ss
			# if count not in site_acc.keys():
			# 	site_acc[count] = acc
	#
	# ss_mapping = {'-':1,'H':2, 'S':3, 'G':4, 'T':5, 'E':6, 'B':7, 'I':8}
	# ss = [0 for i in range(len(sites))]  # 有的位点在dssp文件中没有 对应信息设为0
	# acc = [0 for i in range(len(sites))]
	#
	# for i in range(len(sites)):
	# 	if (i+1) in site_ss.keys():
	# 		ss[i] = ss_mapping[site_ss[i+1]]
	#
	#
	# for i in range(len(sites)):
	# 	if (i+1) in site_acc.keys():
	# 		acc[i] = site_acc[i+1]
	# print(len(acc))

	new_file(pdb_name,chain,'secondary_structure', ss)
	add_content(pdb_name,chain,'ACC', acc)


if __name__ == '__main__':

	pdb_path = r'03del_low_pldd'
	error=[]
	infor = pd.read_excel(r'06all_infor.xlsx')
	for i in infor['ACC_ID'].unique():
		pdb_file_name = i +'.pdb'
		pdb_chain = 'A'
		pdb_file = os.path.join(pdb_path, pdb_file_name)
		get_second_structure_acc(pdb_file, pdb_chain)
		# break

		# print(pdb_file[-10:-6])
		# run_dssp(pdb_file)

# 		except:
# 			error.append(row['pdb'])
# 	print(error)
# 	print(len(error))
# 	pdb_file =r'D:\ruanjian\wendang\pycharm1\PTM_FUNCRION_CLEAN\feature\final_dele_hetatm\1NHL_A.pdb'
# 	# # binding_file = '/home/zjliang/users/zgy/capsule-network/calc_feature/final_pdb/1a4m_pocket.pdb'
# 	pdb_chain = 'A'
# 	# pdb_site, binding_site = get_pdb_binding_site(pdb_file, binding_file, chain)
# 	get_second_structure_acc(pdb_file, pdb_chain)








