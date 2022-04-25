import streamlit as st
import pandas as pd
import math 
#_________________________
import py3Dmol
from stmol import showmol
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
#_________________________
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D

from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem import rdMolDescriptors
#_________________________
#Inicio#
st.title ("FENÓMENOS CUÁNTICOS")
st.subheader("Erick López Saldviar 348916")





seleccion = st.selectbox("Seleccione una opción: ", ["Visualizacion molecular", "Reactividad", "Otros parametros"])

#__________________________________________________________________________________________
#Calculo reactividad#
if seleccion == "Reactividad":
  with st.form(key='calc_react'):
    st.write("Bienvenido. Este programa te calculara parámetros de reactividad")
    st.subheader('Hartress')
    ht = st.number_input("Hartress: ")
    ht0 = st.number_input("Hartress 0: ")
    ht1p = st.number_input("Hartress +1: ")
    ht1m = st.number_input("Hartress -1: ")
    st.subheader('Nucleofilicidad')
    homo = st.number_input("HOMO:")
    lumo = st.number_input("LUMO: ")
    homo_o = st.number_input("HOMO (O):")
    lumo_o = st.number_input("LUMO (V): ")

    calcular = st.form_submit_button('Calcular')
 
 
  if calcular:
    eV0 = ht0 * 27.2116
    eV1p = ht1p * 27.2116
    eV1m = ht1m * 27.2116
  
    kcal = str(ht*627.5)

    A = eV0-eV1m
    I = eV1p - eV0

    n = (I-A)/2
    u = (I+A)/2
    w = pow(u,2)/(2*n)


    col1, col2, col3, col4 = st.columns(4)
    col1.metric(label="H 0", value=ht0)
    col2.metric(label="H -1", value=ht1m)
    col3.metric(label="H +1", value=ht1p)
    col4.metric(label="kCal/mol", value=kcal)

    col1, col2, col3 = st.columns(3)
    col1.metric(label="eV 0", value=str(eV0))
    col2.metric(label="eV -1", value=str(eV1p))
    col3.metric(label="eV +1", value=str(eV1m))
  
    st.subheader("Aproximación de energias")

    col1, col2 = st.columns(2)
    col1.metric(label="Afinidad electrónica", value=str(A))
    col2.metric(label="Potencial de ionización", value=str(I))

    col1, col2, col3 = st.columns(3)
    col1.metric(label="Dureza", value=str(n))
    col2.metric(label="Electronegatividad", value=str(u))
    col3.metric(label="Electrofilicidad", value=str(w))

    st.subheader("Aproximación Orbital")

    A_orb = -27*lumo
    I_orb = -27*homo

    n_orb = (I_orb - A_orb) / 2
    u_orb = (I_orb + A_orb) / 2
    w_orb = pow(u, 2) / (2 * n)
  
    col1, col2 = st.columns(2)
    col1.metric(label="Afinidad electrónica", value=str(A_orb))
    col2.metric(label="Potencial de ionización", value=str(I_orb))

    col1, col2, col3 = st.columns(3)
    col1.metric(label="Dureza", value=str(n_orb))
    col2.metric(label="Electronegatividad", value=str(u_orb))
    col3.metric(label="Electrofilicidad", value=str(w_orb))

    gap = lumo_o - homo_o
    gap_eV = gap*27.2116
  
    col1, col2 = st.columns(2)
    col1.metric(label="GAP", value=str(gap))
    col2.metric(label="GAP (eV", value=str(gap_eV))
#__________________________________________________________________________________________
#Visualización molecular#
if seleccion == "Visualizacion molecular":
    st.title('VISUALIZACIÓN MOLECUALR')
    st.write("Bienvenido. Aquí podrás ver la molecula en su forma tridimensional")
    
    seleccion_molecula = st.selectbox("Seleccione una opción: ", ["SMILES", "Subir un archivo"])
    if seleccion_molecula == "SMILES":
      compound_smiles=st.text_input('SMILES please','COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21')

    if seleccion_molecula == "Subir un archivo":
      uploaded_files = st.sidebar.file_uploader("Choose xyz files", accept_multiple_files=True)

    for uploaded_file in uploaded_files:
      xyz = uploaded_file.getvalue().decode("utf-8")
      render_mol(xyz)

    def render_mol(xyz):
      xyzview = py3Dmol.view(width=400,height=400)
      xyzview.addModel(xyz,'xyz')
      xyzview.setStyle({'stick':{}})
      xyzview.setBackgroundColor('white')#('0xeeeeee')
      xyzview.zoomTo()
      showmol(xyzview, height = 500,width=800)

    def makeblock(smi):
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        mblock = Chem.MolToMolBlock(mol)
        return mblock

    def render_mol(xyz):
        xyzview = py3Dmol.view()#(width=1000,height=1000)
        xyzview.addModel(xyz,'mol')
        xyzview.setStyle({'stick':{}})
        xyzview.setBackgroundColor('white')
        xyzview.zoomTo()
        showmol(xyzview,height=500,width=500)

    
    blk=makeblock(compound_smiles)
    render_mol(blk)
#__________________________________________________________________________________________
#Otros parámetros#
if seleccion == "Otros parametros":
    st.write("Bienvenido. Aqui encontraras LogP y otras cosas")
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
