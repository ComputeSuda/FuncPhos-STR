import os
def run_dssp(pdb_file):
	pdb_name = pdb_file[-10:-4]
	output_path = 'dssp_out/'
	if os.path.exists(output_path) == False:
		os.makedirs(output_path)
	output_file = output_path + pdb_name + '_dssp.out'
	os.system('dssp -i ' + pdb_file + ' -o ' + output_file)  # 在pdb文件上运行dssp

if __name__ == '__main__':
	path='01acc'
	file=os.listdir(path)
	for i in file:
		name=os.path.join(path,i)
		run_dssp(name)