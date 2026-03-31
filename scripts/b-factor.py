#!/usr/bin/env python3

import os
import pandas as pd 
from biopandas.pdb import PandasPdb 
from biopandas.mmcif import PandasMmcif # biopandas cannot write cif files, try exploring CIF FileWriter or Gemmi????
import gemmi # cif parser, use to save cif file
import sys

print("Python executable:", sys.executable)

# test: 2F1W

def get_structure(pdb_id, csp_file): 
    'Load user inputted structure ID and NMR data file'

    try:
        structure = PandasPdb().fetch_pdb(pdb_id)
        if structure is None:
            print("No valid structure loaded.")
            return None, None
        print("Loaded structure from PDB ID:", pdb_id)
    except Exception as e:
        print("Error fetching PDB ID:", e)
        return None, None
    
    csp_df = read_csp_file(csp_file)
    return structure, csp_df

def detect_filetype(filepath): 
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
    
def read_csp_file(csp_file): 
    try:
        df = pd.read_excel(csp_file, sheet_name=1)
        print("Loaded CSP data from Excel.")
    except Exception as e:
        print("Error loading CSP file:", e)
        return None, None
    csp_df = df.rename(columns={"Unnamed: 0" : "Residue", "1to4" : "CSP"}) # rename column headers
    csp_df.columns = ['Residue', 'CSP']
    csp_dict = dict(zip(csp_df['Residue'], csp_df['CSP'])) 
    return csp_dict
    
def map_file(filepath, csp_df):  
    '''If using user-provided structure file instead of PDB ID'''

    structure, filetype = detect_filetype(filepath)
    if structure is None:
        print("Error: could not load structure file.")
        return None

    atom_df=structure.df['ATOM'].copy() # get data with atom key
    print(structure.df['ATOM'].columns) # check key names in file 
    
    if filetype == "PDB":
        res_num_col = 'residue_number'
        res_name_col = 'residue_name'
        b_factor_col = 'b_factor'  
        # print(atom_df[['residue_number', 'residue_name', 'b_factor']].head(10)) # check before refactoring

    if filetype == "CIF":
        
        res_num_col = 'auth_seq_id'
        res_name_col = 'auth_comp_id'
        b_factor_col = 'B_iso_or_equiv'

        # change name of keys to be same as PDB file
        atom_df['residue_number'] = atom_df[res_num_col].astype(int) # try
        atom_df['residue_name'] = atom_df[res_name_col]
        atom_df['b_factor'] = atom_df[b_factor_col].astype(float) 

    return atom_df, structure, filetype

def substitute_b_factor_using_id(structure_id, csp_df, save_dir=None):  
    '''Substitutes b-factor column with chemical shift perturbation data using a the PDB ID; will output PDB file'''
    structure, csp_df = get_structure(structure_id, csp_df)

    if structure is None:
        return None, None
    
    atom_df=structure.df['ATOM'].copy() # get data with atom key

    atom_df['b_factor'] = atom_df['residue_number'].map(csp_df).fillna(0.0) # map values to residues


    # put results back into data frame
    structure.df['ATOM'] = atom_df 

    # define save directory
    if save_dir is None:
        save_dir = os.getcwd()
    out_filename = f"{structure_id}_b-factor.pdb"
    out_path = os.path.join(save_dir, out_filename)
    
    structure.to_pdb(path=out_path, records=['ATOM'], gz=False)
    print(f"Saved PDB with substituted B-factors to: {out_path}")
    return out_path


def substitute_b_factor_using_file(filepath, csp_df): 
    '''Substitutes b-factor column with chemical shift perturbation data using a user provided protein file; will output PDB or CIF file'''
    structure, filetype = detect_filetype(filepath)
    if structure is None:
        print("No valid structure loaded.")
        return None
    
    atom_df = structure.df['ATOM'].copy()

    csp_dict = dict(zip(csp_df['Residue'], csp_df['CSP'])) # convert csp pandas df to file to prepare for merge
    atom_df['b_factor'] = atom_df['residue_number'].map(csp_dict).fillna(0.0)

    if filetype == "PDB":
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

def prompt_user(): 
    '''Prompts user to input csp data, protein id or filetype, and save directory'''
    csp_file = input("Input the file path to your data file: ") 
    csp_df = read_csp_file(csp_file)
    
    if csp_df is None:
        print("Error reading CSP file. Exiting.")
        return
    

    save_dir = input("Input directory to save output (or press enter for current dir): ") # Pass this through later in function so that it subs None if they choose one
    save_dir = save_dir if save_dir.strip() else None

    user_choice = input("Would you like to upload a structure file? (Type 'Y' or 'N'): ") 

    if user_choice.upper() == "Y": 
        structure_file = input("Please provide a path to either a PDB or CIF file: ")
        atom_df, structure, filetype = map_file(structure_file, csp_df)
        out_path = substitute_b_factor_using_file(structure_file, csp_df)
        return out_path, f'B-factor substitution complete. Your new file can be found at {out_path}.'
    
    if user_choice.upper() == "N": 
        pdb_id = input("Input a PDB ID: ")   
        structure, out_path = substitute_b_factor_using_id(pdb_id, csp_df)
        return out_path, f'B-factor substitution complete. Your new file can be found at {out_path}.'

    return user_choice, csp_file, save_dir

def main(): 
    user_choice, csp_file, save_dir = prompt_user()

if __name__ == "__main__":
    main()


# TODO dubug keys in csp_dict at 119, KeyError: 'Residue' - likely that csp_file/csp_df/csp_dict is mixed up. Check filetypes before getting to line 117 to see what the status of the csp object is. 