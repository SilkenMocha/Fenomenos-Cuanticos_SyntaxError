import streamlit as st
import py3Dmol
from stmol import showmol
import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem import rdMolDescriptors

compound_smiles=st.text_input('SMILES please','COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21')

st.title(compound_smiles)
m = Chem.MolFromSmiles(compound_smiles)
tpsa = Descriptors.TPSA(m)
logP = Descriptors.MolLogP(m)
st.write(tpsa)
st.write(logP)


mol = Chem.MolFromSmiles(compound_smiles)
# Gasteiger partial charges
AllChem.ComputeGasteigerCharges(mol)
contribs = [mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge') for i in range(mol.GetNumAtoms())]
fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, contribs, colorMap='jet', contourLines=10)
st.pyplot(fig)

# Crippen contributions to logP
contribs = rdMolDescriptors._CalcCrippenContribs(mol)
fig2 = SimilarityMaps.GetSimilarityMapFromWeights(mol,[x for x,y in contribs], colorMap='jet', contourLines=10)
st.pyplot(fig2)
