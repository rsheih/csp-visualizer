#!/usr/bin/env python3

# SCRIPT USED WHEN USER HAS PDB OR CIF FILE

import os
import pandas as pd 
from biopandas.pdb import PandasPdb
from biopandas.mmcif import PandasMmcif # biopandas cannot write cif files, try exploring CIF FileWriter or Gemmi????
import gemmi # cif parser, use to save cif file
import sys
print(sys.executable)

#######################
# OPEN STRUCTURE FILE #
#######################

def specify_file(filepath): 
    '''Detects file type and returns as variable'''
    ext = os.path.splitext(filepath)[1].lower()
    if ext == ".pdb":
        ppdb=PandasPdb().read_pdb(filepath)
        return ppdb, "PDB"
    elif ext in [".cif", ".mmcif"]:
        ppdb=PandasMmcif().read_mmcif(filepath)
        return ppdb, "CIF"
    else:
        return None, None

####################
# CREATE NEW FILE #
###################

def substitute_b_factor(filepath, csp_df):  
    '''Substitutes b-factor column with chemical shift perturbation data; will output PDB file'''
    structure, filetype = specify_file(filepath) 
    # csp_df = pd.read_excel(csp_file) # excel file with chem shift info used to substitute b factor values

    if structure is None:
        print("No valid structure loaded. Exiting.")
        return None
    
    # rename column headers
    csp_df.columns = ['Residue', 'CSP']
    atom_df=structure.df['ATOM'].copy() # get data with atom key
    print(structure.df['ATOM'].columns) # check key names in file 
    
    # map key names in filetype for either PDB or CIF
    if filetype == "PDB":
        res_num_col = 'residue_number'
        res_name_col = 'residue_name'
        b_factor_col = 'b_factor'  
        print(atom_df[['residue_number', 'residue_name', 'b_factor']].head(10)) # check before refactoring
        csp_dict = dict(zip(csp_df['Residue'], csp_df['CSP'])) # convert csp pandas df to file to prepare for merge
        atom_df['b_factor'] = atom_df['residue_number'].map(csp_dict).fillna(0.0) # map values to residues

        # check after merge
        check_df = atom_df[['residue_number', 'residue_name', 'b_factor']] 
        print(check_df.head(20))

        # put results back into data frame
        structure.df['ATOM'] = atom_df 

        # define output file path
        dirname, filename = os.path.split(filepath) 
        name = os.path.splitext(filename)
        out_path = os.path.join(dirname, f"TEST_{name}_b-factor.pdb") 
        print(f"Saved to {out_path}")
        return structure.to_pdb(path=out_path, records=['ATOM'], gz=False) 
    
    if filetype == "CIF":
        
        # set keys
        
        res_num_col = 'auth_seq_id'
        res_name_col = 'auth_comp_id'
        b_factor_col = 'B_iso_or_equiv'

        # change name of keys to be same as PDB file
        atom_df['residue_number'] = atom_df[res_num_col].astype(int) # try
        atom_df['residue_name'] = atom_df[res_name_col]
        atom_df['b_factor'] = atom_df[b_factor_col].astype(float)

        # convert excel df to dict to prepare for merge
        csp_dict = dict(zip(csp_df['Residue'], csp_df['CSP'])) 

        # set new key ['new_b']
        atom_df['new_b'] = atom_df['residue_number'].map(csp_dict).fillna(0.0)

        print(atom_df[['residue_number', 'new_b']].head(10))  # debug

        print("inspecting object")

        # read cif file using gemmi 
        doc = gemmi.cif.read_file(filepath)
        block = doc.sole_block()

        # configure data structure 
        loop = None
        for item in block:
            if item.loop is not None and '_atom_site.id' in item.loop.tags:
                loop = item.loop
                break

        if loop is None:
            print("ERROR: atom_site loop not found")
            return
        
        # find b factor index in cif file
        b_idx = None
        for i, tag in enumerate(loop.tags):
            if 'B_iso_or_equiv' in tag:
                b_idx = i
                print("Found B-factor column:", tag) 
                break
        
        # set value for b factor to be equal to new-z
        for i in range(len(atom_df)):
            loop[i, b_idx] = f"{atom_df.iloc[i]['new_b']:.2f}"

        for tag in atom_df.columns:
            if 'B_iso_or_equiv' in tag:
                atom_df[tag] = atom_df['residue_number'].map(csp_dict).fillna(0.0)
                print("Replaced:", tag)
       
       # save
        dirname, filename = os.path.split(filepath)
        name = os.path.splitext(filename)[0]
        out_path = os.path.join(dirname, f"TEST_{name}_b-factor.cif")

        print(f"Saved CIF to {out_path}")
        doc.write_file(out_path)

        return out_path

def main():  
    '''Substitutes b factor column in user input PDB/CIF file 
    with chemical shift data from XLSX file and outputs a PDB file'''
    filepath = input("Enter file path: ") 

    # csp_df = pd.read_excel("/Users/rebekahsheih/Library/CloudStorage/OneDrive-UConnHealthCenter/rotation_3_bezsonova/usp7_files/USP7-SCML2-NMR-titration.xlsx",  sheet_name=1)  
    csp_df=pd.read_excel("/Users/rebekahsheih/projects/bezsonova_lab/csp-visualizer/usp7_files/USP7-SCML2-NMR-titration.xlsx", sheet_name=1)
    # TODO (Use text file with two columns (residue and chem shift) instead of excel file)
    print('loaded excel')
    # csp_df=input("Enter a .txt file: ")
    new_structure = substitute_b_factor(filepath, csp_df)
    return new_structure

if __name__ == "__main__": 
    main()

# TEST 
# filepath=('/Users/rebekahsheih/projects/bezsonova_lab/usp7_files/2F1W.pdb')

# TODO 
# Add option for them to choose where path is saved 
# Add option for inputting .txt file  
# import dependencies so that they have the right packages installed