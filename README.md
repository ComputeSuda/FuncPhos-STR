
# FuncPhos-STR
FuncPhos-STR:The model, which consists of three feature encoding sub-networks (StrNet, SeqNet and SpNet) and a heterogeneous feature combination sub-network CoNet, comprehensively uses protein structure, protein sequence and PPI network information to predict functional phosphosites.

# System requirement
FuncPhos-SEQ is develpoed under Linux environment with:
* Python (3.10.4):
    - keras==2.8.0
    - networkx==2.6.3
    - scipy==1.7.3
    - scikit-learn==0.24.2
    - numpy==1.23.3
    - tensorflow==2.8.0
    - biopython==1.78
    - prody==2.0
* You can install the dependent packages by the following commands:
    - pip install python==3.10.4
    - pip install numpy==1.23.3
    - pip install keras==2.8.0
    - pip install tensorflow==2.8.0
# Dataset
We provide phosphosite data, collected from five databases - PSP, EPSD, PLMD, IPTMNet and PTMD - detailing information on phosphosites and their regulation of molecular functions, biological processes and intermolecular interactions.
