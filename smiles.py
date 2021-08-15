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

################ Sidebar ####################
with st.sidebar.expander('Rule One'):
  st.markdown('''
## Atoms
|If |then |
|----|----|
| Non-aromatic atoms |Uper case letters |
| Aromatic atoms |lower case letters |
|Atomic symbols has more than one letter | The second is lower case |
## Bonds
| Bond type| Bond symbol |
|---|---|
|Simple | - |
|Double|=|
|Triple|#|
|Aromatic|*|
| Disconnected structures|. |

### Example:
 CC   👉 There is a non-aromatic carbon attached to another non-aromatic carbon by a single bond.

🛑 A bond between two lower case atom symbols is *aromatic*.
''')

with st.sidebar.expander('Rule Two'):
  st.markdown('''
  ## Simple chains
  * Structures are hydrogen suppresed (Molecules represented without hydrogens)
  * If enough bonds are not identified by the user, the system will assume that connections
  are satisfied by hidrogens.
  * The user can explicitly identify hydrogen bonds, but if so the interpreter will assume that all of them are fully identified.
  Note:
  
  *Because SMILES allows entry of all elements in the periodic table, 
  and also utilizes hydrogen suppression, the user should be aware of chemicals with two letters 
  that could be misinterpreted by the computer. For example, 'Sc' could be interpreted as a **sulfur**
  atom connected to an aromatic **carbon** by a single bond, or it could be the symbol for **scandium**. 
  The SMILES interpreter gives priority to the interpretation of a single bond connecting a sulfur atom and an aromatic carbon. 
  To identify scandium the user should enter [Sc]*.
  ''')

with st.sidebar.expander('Rule Three'):
  st.markdown('''
  ## Branches
  * A branch from a chain is specified by placing the SMILES symbol(s) for the branch between parenthesis. 
  * The string in parentheses is placed directly after the symbol for the atom to which it is connected. 
  * If it is connected by a double or triple bond, the bond symbol immediately follows the left parenthesis.
  ''')

with st.sidebar.expander('Rule Four'):
  st.markdown('''
  In progress ...
  ''')

with st.sidebar.expander('Rule Five'):
  st.markdown('''
  In progress ...
  ''')
st.sidebar.markdown('Author: José Manuel Nápoles ([@napoles3d](https://twitter.com/napoles3D))')
st.sidebar.write('Info about SMILES: https://archive.epa.gov/med/med_archive_03/web/html/smiles.html')
