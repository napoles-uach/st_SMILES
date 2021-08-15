import streamlit as st
import streamlit.components.v1 as components
import py3Dmol
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

st.title('SMILES  + RDKit + Py3DMOL :smiley:')
def show(smi, style='stick'):
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    mblock = Chem.MolToMolBlock(mol)

    view = py3Dmol.view(width=400, height=400)
    view.addModel(mblock, 'mol')
    view.setStyle({style:{}})
    view.zoomTo()
    view.show()
    view.render()
    t =view.js()
    f = open('viz.html', 'w')
    f.write(t.startjs)
    f.write(t.endjs)
    f.close()

compound_smiles=st.text_input('SMILES please','CC')
m = Chem.MolFromSmiles(compound_smiles)

Draw.MolToFile(m,'mol.png')


show(compound_smiles)
HtmlFile = open("viz.html", 'r', encoding='utf-8')
source_code = HtmlFile.read() 
c1,c2=st.columns(2)
with c1:
  st.write('Molecule :coffee:')
  st.image('mol.png')
with c2:
  components.html(source_code, height = 400,width=400)
