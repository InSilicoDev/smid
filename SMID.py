#!/usr/bin/python3
import sys
import os
import rdkit 
import pandas as pd
from rdkit.Chem import PandasTools
from rdkit import RDConfig
from rdkit import Chem
from rdkit.Chem import AllChem
def SMID(filename):
    multi_sdf = Chem.SDMolSupplier(filename, removeHs=False)
    mxdi=[]
    for mol in multi_sdf:
        mxdi.append(float(AllChem.Get3DDistanceMatrix(mol).max()))
    sdfFile_pd = (filename)
    mol_df = PandasTools.LoadSDF(sdfFile_pd, molColName='Molecule')
    mol_df["mxdi"] = mxdi
    rdkit.Chem.PandasTools.WriteSDF(mol_df, f'output_{filename}', molColName='Molecule', idName=None, properties=list(mol_df.columns))
filename = '' 
if len(sys.argv) > 1:
    filename = sys.argv[1]
    print("Operation Successful")
else:
    print("Unsupported Arguments")
if __name__ == '__main__':
    SMID(filename)