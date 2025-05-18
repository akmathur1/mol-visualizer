import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import py3Dmol
from io import BytesIO
from PIL import Image

st.set_page_config(page_title="Molecule Structure Visualizer by Arjun Mathur", layout="wide")

st.title("🔬 Molecule Structure Visualizer")

# SMILES Input
smiles = st.text_input("Enter a SMILES string (e.g. CC(=O)OC1=CC=CC=C1C(=O)O):", "CCO")

try:
    mol = Chem.MolFromSmiles(smiles)
    mol_3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_3d)
    AllChem.UFFOptimizeMolecule(mol_3d)
except Exception as e:
    st.error(f"❌ Invalid SMILES string or 3D structure generation failed:\n{e}")
    mol = None

if mol:
    # 2D Rendering
    st.subheader("2D Structure")
    img = Draw.MolToImage(mol, size=(300, 300))
    buf = BytesIO()
    img.save(buf, format="PNG")
    st.image(Image.open(buf), caption="RDKit 2D Rendering")

    # 3D Viewer
    st.subheader("3D Interactive Viewer")
    mb = Chem.MolToMolBlock(mol_3d)
    viewer = py3Dmol.view(width=500, height=400)
    viewer.addModel(mb, "mol")
    viewer.setStyle({"stick": {}})
    viewer.setBackgroundColor("white")
    viewer.zoomTo()

    # ✅ Streamlit-compatible render
    st.components.v1.html(viewer._make_html(), height=400)
