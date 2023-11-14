'''activity'''

import os
from prody import *
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
# import structuremap.utils
# structuremap.utils.set_logger()
from processing import download_alphafold_cif, download_alphafold_pae, format_alphafold_data, annotate_accessibility, get_smooth_score
import re
# import plotly.express as px
# import tqdm
# import tempfile
standard_amino_acids = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
						'GLY': 'G', 'GLN': 'Q', 'GLU': 'E', 'HIS': 'H', 'ILE': 'I',
						'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PRO': 'P', 'PHE': 'F',
						'SER': 'S', 'THR': 'T', 'TYR': 'Y', 'TRP': 'W', 'VAL': 'V'}

def new_file(pdb_name, name, data):
    dataframe = pd.DataFrame({name: data})
    csv_name = pdb_name + '_feature.csv'
    dataframe.to_csv(csv_name, index=False)
def add_content(pdb_name, name, data):
    csv_name = pdb_name + '_feature.csv'
    file = pd.read_csv(csv_name)
    # file.insert(index_num,name,data)
    file[name] = data
    file.to_csv(csv_name, index=False)
def process_pdb(pdb_file,pdb_del_plddt_file):
    sequence = ''
    site_sequence_no = []
    site_no = []
    with open(pdb_del_plddt_file,'w') as f:
        with open(pdb_file, 'r') as f1:
            for line in f1:
                if line[:4]!='ATOM':
                    f.write(line)
                else:
                    sitenum = int(line[22:28].replace(' ', ''))
                    residue = line[17:20].replace(' ', '')
                    if sitenum not in site_sequence_no:
                        site_sequence_no.append(sitenum)
                        sequence += standard_amino_acids[residue]
                    line_list=line.strip().split()
                    if float(line_list[-2])>50:
                        f.write(line)
                        if sitenum not in site_no:
                            site_no.append(sitenum)

                    else:
                        continue
    print(len(sequence))
    print(len(site_no))
    new_file(pdb_file[:-4], 'position', site_no)
    return sequence,site_no

    
def change_value(cc):
    # 
    for i in range(cc.shape[0]):
        for j in range(cc.shape[1]):
            if cc[i, j] < 0:
                cc[i, j] = -cc[i, j]
def get_anm_cc(pdb_file):  # binding_site,
    pdb_name = pdb_file[:6]
    structure = parsePDB(pdb_file)
    choose = 'not ion and name CA and chain ' + 'A'
    calphas = structure.select(choose)
    n_modes_n = None  # 
    anm = ANM(pdb_name)
    anm.buildHessian(calphas, cutoff=15.0, gamma=1)
    anm.calcModes(n_modes=n_modes_n)
    print('cal anm_stiffness...')
    anm_sti = calcMechStiff(anm, calphas)  # (349, 349)
    anm_stiffness = np.mean(anm_sti, axis=1)  #
    add_content(pdb_name, 'anm_stiffness', anm_stiffness)

    print('cal anm_cc...')
    n_modes = anm._n_modes
    anm_greater_60_per_cc = calcCrossCorr(anm[int(n_modes * 0.6):])
    change_value(anm_greater_60_per_cc)
    anm_greater_60_per_cc = np.mean(anm_greater_60_per_cc, axis=1)
    add_content(pdb_name, 'anm_greater_60_per_cc', anm_greater_60_per_cc)

    print('cal gnm_eigenvectors...')
    gnm = GNM(pdb_name)
    gnm.buildKirchhoff(calphas, cutoff=10.0, gamma=1)
    gnm.calcModes(n_modes=n_modes_n)
    gnm_top20_eigenvectors = gnm[:20].getEigvecs()

    for i in range(gnm_top20_eigenvectors.shape[0]):
        for j in range(gnm_top20_eigenvectors.shape[1]):
            if gnm_top20_eigenvectors[i, j] < 0:
                gnm_top20_eigenvectors[i, j] = -gnm_top20_eigenvectors[i, j]
    gnm_top20_eigenvectors = np.mean(gnm_top20_eigenvectors, axis=1)
    add_content(pdb_name,'gnm_top20_eigenvectors', gnm_top20_eigenvectors)


def get_nacen_feature(r_path,R_path,feature_path):
    os.system(r_path+' '+R_path)
    pdb_name = pdb_file[:-4]
    nacen_feature = pd.read_csv(feature_path+pdb_name+'_nacen.csv')
    # print(len(nacen_feature))
    Unweighted_Closeness=[]
    Node_weighted_degree=[]
    Unweighted_Degree=[]
    Node_weighted_Closeness=[]
    for i in range(len(nacen_feature)):
        site_feature = list(nacen_feature.iloc[i, :])[0].split('\t')
        # print(site_feature)
        # break
        Unweighted_Closeness.append(site_feature[6])
        Unweighted_Degree.append(site_feature[4])
        Node_weighted_degree.append(site_feature[7])
        Node_weighted_Closeness.append(site_feature[9])

    # print(Unweighted_Closeness)
    add_content(pdb_name, 'Unweighted_Closeness', Unweighted_Closeness)
    add_content(pdb_name, 'Node_weighted_degree', Node_weighted_degree)
    add_content(pdb_name, 'Unweighted_Degree', Unweighted_Degree)
    add_content(pdb_name, 'Node_weighted_Closeness', Node_weighted_Closeness)
    # site_feature = list(nacen_feature.iloc[0:, :])
    # print(nacen_feature.loc[0])
    # print(site_feature)
def get_structuremap_feature(pdb_name,site_no):
    test_proteins = [pdb_name]
    nAA_24_180_pae_smooth10=[]
    # secondery_structure=[]
    # structure_map={'unstructured':0,'HELX_LH_PP_P':1,'HELX_RH_AL_P':2,'HELX_RH_3T_P':3,'STRN':4,'BEND':5,'TURN_TY1_P':6}
    cif_dir = os.path.join('alpha_pdb_file', 'tutorial_cif')
    pae_dir = os.path.join('alpha_pdb_file', 'tutorial_pae')
    # '''download'''
    valid_proteins_cif, invalid_proteins_cif, existing_proteins_cif = download_alphafold_cif(
        proteins=test_proteins,
        out_folder=cif_dir)
    valid_proteins_pae, invalid_proteins_pae, existing_proteins_pae = download_alphafold_pae(
        proteins=test_proteins,
        out_folder=pae_dir,
        )
    alphafold_annotation = format_alphafold_data(
        directory=cif_dir,
        protein_ids=test_proteins)
    '''#ppse'''
    full_sphere_exposure = annotate_accessibility(
        df=alphafold_annotation,
        max_dist=24,
        max_angle=180,
        error_dir=pae_dir)
    alphafold_accessibility = alphafold_annotation.merge(
        full_sphere_exposure, how='left', on=['protein_id', 'AA', 'position'])
    alphafold_accessibility_smooth = get_smooth_score(
        alphafold_accessibility,
        np.array(['nAA_24_180_pae']),
        [10])
    for index,row in alphafold_accessibility_smooth.iterrows():
        # print(row['position'])
        if row['position'] in site_no:
            nAA_24_180_pae_smooth10.append(row['nAA_24_180_pae_smooth10'])
    print(nAA_24_180_pae_smooth10)
    print(len(nAA_24_180_pae_smooth10))
    add_content(pdb_name, 'nAA_24_180_pae_smooth10', nAA_24_180_pae_smooth10)

def get_iij_feature(r_path,R_path,feature_path):
    iij_feature = feature_path + pdb_name + '_iij.txt'
    if os.path.exists(iij_feature) == False:
        os.system(r_path + ' ' + R_path)
    closeness = []
    page_rank = []
    i = 0

    with open(iij_feature, 'r') as f:
        for line in f:
            if line.strip().replace('"', '') == 'closeness':
                i = 2
                continue
            if line.strip().replace('"', '') == 'page_rank':
                i = 6
                continue
            if i == 2:
                closeness.append(line.strip().split(' ')[-1])
            if i == 6:
                page_rank.append(line.strip().split(' ')[-1])
    # print(len(closeness))
    add_content(pdb_name, 'closeness', closeness)
    add_content(pdb_name, 'page_rank', page_rank)
def get_second_acc_feature(r_path,R_path,feature_path):
    # os.system(r_path + ' ' + R_path)
    pdb_name = pdb_file[:-4]
    dssp_feature = pd.read_csv(feature_path + pdb_name + '_dssp.csv')
    # print(len(nacen_feature))
    secondery_structure = []
    ACC = []
    secondery_map={'X':0,'B':1,'E':2,'G':3,'H':4,'I':5,'T':7,'S':8}
    for i in range(len(dssp_feature)):
        site_feature = list(dssp_feature.iloc[i, :])[0].split()
        # print(site_feature)
        secondery_structure.append(secondery_map[site_feature[1]])
        ACC.append(site_feature[2])
        # print(secondery_structure)
    # print(Unweighted_Closeness)
    add_content(pdb_name, 'secondery_structure', secondery_structure)
    add_content(pdb_name, 'ACC', ACC)

def get_site_feature(pdb_name,site_list,model_name):
    csv_name = pdb_name + '_feature.csv'
    file=pd.read_csv(csv_name)
    str_feature=[]
    activity_feature=['anm_greater_60_per_cc', 'Unweighted_Closeness','gnm_top20_eigenvectors','Node_weighted_degree','anm_stiffness','nAA_24_180_pae_smooth10']
    ppi_feature=['anm_greater_60_per_cc','gnm_top20_eigenvectors','page_rank','Unweighted_Degree','Node_weighted_Closeness','ACC']
    function_feature=['anm_greater_60_per_cc','closeness','gnm_top20_eigenvectors','page_rank','Unweighted_Degree','secondery_structure']
    index_list=[]
    if model_name=='Activity':
        for f in activity_feature:
            index_list.append(file.columns.get_loc(f))
    elif model_name=='Function':
        for f in function_feature:
            index_list.append(file.columns.get_loc(f))
    else:
        for f in ppi_feature:
            index_list.append(file.columns.get_loc(f))
    print(index_list)
    for index,row in file.iterrows():
        if row['position'] in site_list:
            str_feature.append(row[index_list])
    str_feature=np.array(str_feature)
    print(str_feature.shape)
    str_file_name=pdb_name+'_'+model_name+'_str.npy'
    np.save(str_file_name,str_feature)

def get_model_feature(model_name,pdb_name,pdb_file,pdb_del_plddt_file,site_list):
    r_path = r'Rscript.exe'
    R_nacen_path = r'NACEN_f.R'
    R_iij_path=r'iij_test.R'
    R_dssp_path=r'dssp.R'
    feature_path = 'feature/'

    if model_name=='Activity':
        # sequence, site_no = process_pdb(pdb_file, pdb_del_plddt_file)
        # get_anm_cc(pdb_del_plddt_file)
        # get_nacen_feature(r_path, R_nacen_path, feature_path)
        # get_structuremap_feature(pdb_name, site_no)
        get_site_feature(pdb_name, site_list,model_name)

    elif model_name=='Function':
        sequence, site_no = process_pdb(pdb_file, pdb_del_plddt_file)
        get_anm_cc(pdb_del_plddt_file)
        get_nacen_feature(r_path, R_nacen_path, feature_path)
        get_iij_feature(r_path,R_iij_path,feature_path)
        get_second_acc_feature(r_path,R_dssp_path,feature_path)
        get_site_feature(pdb_name, site_list, model_name)
    else:
        sequence, site_no = process_pdb(pdb_file, pdb_del_plddt_file)
        get_anm_cc(pdb_del_plddt_file)
        get_nacen_feature(r_path, R_nacen_path, feature_path)
        get_iij_feature(r_path, R_iij_path, feature_path)
        get_second_acc_feature(r_path, R_dssp_path, feature_path)
        get_site_feature(pdb_name, site_list, model_name)


'''ppi:anm_greater_60_per_cc;gnm_top20_eigenvectors;page_rank;
Unweighted Degree;Node-weighted Closeness;secondery_structure'''

'''activity:anm_greater_60_per_cc; Unweighted Closeness; gnm_top20_eigenvectors;
Node-weighted degree;anm_stiffness;nAA_24_180_pae_smooth10'''

'''function:anm_greater_60_per_cc;closeness;gnm_top20_eigenvectors;page_rank;
Unweighted Degree;ACC'''

if __name__ == '__main__':
    pdb_name='A1X283'
    site_list=[3,4,5]
    pdb_file='A1X283.pdb'
    model_name='Activity'
    pdb_del_plddt_file = pdb_file[:-4] + '_new.pdb'
    # get_model_feature(model_name,pdb_name,pdb_file,pdb_del_plddt_file,site_list)
    file=np.load('A1X283_ppi_str.npy')
    print(file.shape)
    print(file[:,:128].shape)


#
# Array.prototype._filter=function(fn){
#     if(typeof fn!='function') return
#     var newArr=[]
#     for(i=0;i<this.length;i++){
#        if(fn(this[i],i,this)){
#         newArr.push(this[i])
#        }
#     }
#     return newArr
