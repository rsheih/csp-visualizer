#!/usr/bin/env python3

import os
import importlib
import subprocess 
import sys

print("Python executable:", sys.executable)

# test: 2F1W

def load_dependencies(): 

    global pd, PandasPdb, PandasMmcif, gemmi 

    import pandas as pd 
    from biopandas.pdb import PandasPdb 
    from biopandas.mmcif import PandasMmcif
    import gemmi # cif parser, use to save cif file

def print_user_instructions():
    print("\n" + "="*70)
    print("CSP VISUALIZER — INPUT REQUIREMENTS")
    print("="*70)
    
    print("\n1. Chemical Shift Perturbation (CSP) File Requirements:")
    print("   • File type must be: .txt")
    print("   • Must contain exactly TWO columns")
    print("   • Each column must have a header")
    print("   • Column 1: Residue number")
    print("   • Column 2: Chemical shift perturbation values")
    
    print("\n   Example format:")
    print("   Residue    CSP")
    print("   65          5.68")
    print("   66          6.44")
    
    print("\n2. Structure File OR PDB ID (next step):")
    print("   After uploading your CSP file, you will be prompted")
    print("   to provide a protein structure file.")
    print("   • Accepted formats: .pdb or .cif")
    print("   • You may also provide a PDB ID instead")
    
    print("\n3. Workflow:")
    print("   CSP file → Structure file → B-factor substitution → Output file")
    
    print("\n" + "="*70 + "\n")

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
    
def normalize_csp(df):
    '''Ensure CSP file is always safe for mapping'''

    df = df.iloc[:, :2].copy()
    df.columns = ["Residue", "CSP"]

    df["Residue"] = pd.to_numeric(df["Residue"], errors="coerce")
    df["CSP"] = pd.to_numeric(df["CSP"], errors="coerce")

    df = df.dropna()

    return df
    
def read_csp_file(csp_file): 
    '''reads user input text file'''
    try:
        ext = os.path.splitext(csp_file)[1].lower()

        if ext in [".xlsx", ".xls"]:
            csp_df = pd.read_excel(csp_file)

        elif ext in [".txt", ".tsv"]:
            csp_df = pd.read_csv(csp_file, sep=r"\s+", header=0)  

        else:
            print("Unsupported file type")
            return None

        print("Loaded CSP data.")

    except Exception as e:
        print("Error loading CSP file:", e)
        return None
    
    # csp_df = csp_df.iloc[:, :2].copy()

    csp_df.columns = ["Residue", "CSP"] # overwrite headers

    return normalize_csp(csp_df)
    
def map_keys(filepath):  
    '''Prepare atom_df and rename columns to standard form'''

    structure, filetype = detect_filetype(filepath)
    if structure is None:
        print("Error: could not load structure file.")
        return None, None, None

    atom_df=structure.df['ATOM'] # get data with atom key

    
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
    try:
        structure = PandasPdb().fetch_pdb(structure_id)
        if structure is None:
            print("No valid structure loaded.")
            return None, None
        print("Loaded structure from PDB ID:", structure_id)
    except Exception as e:
        print("Error fetching PDB ID:", e)
        return None

    # dirname, filename = os.path.split(structure_id)
    name = structure_id.strip().upper()
    if save_dir is None:
        save_dir = os.getcwd()

    atom_df=structure.df['ATOM'].copy() # get data with atom key

    print(atom_df[['residue_number', 'residue_name', 'b_factor']].head(10)) # check before refactoring
    csp_dict = dict(zip(csp_df['Residue'], csp_df['CSP'])) # convert csp pandas df to file to prepare for merge
    atom_df['b_factor'] = atom_df['residue_number'].map(csp_dict).fillna(0.0) # map values to residues

    # check after merge
    check_df = atom_df[['residue_number', 'residue_name', 'b_factor']] 
    print(check_df.head(20))

                
    atom_df['b_factor'] = atom_df['residue_number'].map(csp_dict).fillna(0.0)
    structure.df['ATOM'] = atom_df
    out_path = os.path.join(save_dir, f"{name}_b-factor.pdb")

    
    structure.to_pdb(path=out_path, records=['ATOM'], gz=False)
    print(f"Saved PDB with substituted B-factors to: {out_path}")
    return out_path

def substitute_b_factor_using_file(filepath, csp_df): 
    '''Substitutes b-factor column with chemical shift perturbation data using a user provided protein file; will output PDB or CIF file'''
    atom_df, structure, filetype = map_keys(filepath) # returns pandas.DataFrame, pandas.DataFrame, str
    if structure is None:
        print("No valid structure loaded.")
        return None

    if filetype == "PDB":

        atom_df=structure.df['ATOM'].copy() # get data with atom key

        print(atom_df[['residue_number', 'residue_name', 'b_factor']].head(10)) # check before refactoring
        csp_dict = dict(zip(csp_df['Residue'], csp_df['CSP'])) # convert csp pandas df to file to prepare for merge
        atom_df['b_factor'] = atom_df['residue_number'].map(csp_dict).fillna(0.0) # map values to residues

        # check after merge
        check_df = atom_df[['residue_number', 'residue_name', 'b_factor']] 
        print(check_df.head(20))

        # for val in atom_df["b_factor"]: 
        #     print(val)
                
        atom_df['b_factor'] = atom_df['residue_number'].map(csp_dict).fillna(0.0)
        structure.df['ATOM'] = atom_df
        out_path = os.path.join(dirname, f"{name}_b-factor.pdb")
        structure.to_pdb(path=out_path, records=['ATOM'], gz=False)
        print(f"Saved PDB to {out_path}")
        return out_path

    
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
        
        print(atom_df[['residue_number', 'residue_name', 'b_factor']].head(10)) # check before refactoring

        # set new key ['new_b']
        atom_df['new_b'] = atom_df['residue_number'].map(csp_dict).fillna(0.0)


        # print("inspecting object")

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

        # check after merge
        check_df = atom_df[['residue_number', 'residue_name', 'new_b']] 
        print(check_df.head(10))
       
       # save
        dirname, filename = os.path.split(filepath)
        name = os.path.splitext(filename)[0]
        out_path = os.path.join(dirname, f"{name}_b-factor.cif")

        print(f"Saved CIF to {out_path}")
        doc.write_file(out_path)

        return out_path
    
def prompt_user(): 
    '''Prompts user to input csp data, protein id or filetype, and save directory'''

    csp_file = input("Input the file pathway to your Chemical Shift Perturbation (CSP) File: ")  
    # csp_file='/Users/rebekahsheih/projects/bezsonova_lab/csp-visualizer/usp7_files/USP7-SCML2-NMR-titration.xlsx'
    csp_df = read_csp_file(csp_file) # returns df
    
    if csp_df is None:
        print("Error reading CSP file. Exiting.")
        return
    

    save_dir = input("Input directory to save output (or press enter for current dir): ") # Pass this through later in function so that it subs None if they choose one
    save_dir = save_dir if save_dir.strip() else None

    user_choice = input("Would you like to upload a structure file? (Type 'Y' or 'N'): ") 

    if user_choice.upper() == "Y": 
        structure_file = input("Input the file pathway to you Structure File (either a PDB or CIF file): ") 
        # structure_file='/Users/rebekahsheih/projects/bezsonova_lab/csp-visualizer/usp7_files/2F1W.pdb'
        # atom_df, structure, filetype = map_file(structure_file, csp_dict)
        # out_path = substitute_b_factor_using_file(structure_file, csp_dict)
        return substitute_b_factor_using_file(structure_file, csp_df)
    
    if user_choice.upper() == "N": 
        pdb_id = input("Input a PDB ID: ")   
        return substitute_b_factor_using_id(pdb_id, csp_df, save_dir)

def main():
    # check_dependencies() 
    load_dependencies()
    print_user_instructions()
    prompt_user()

if __name__ == "__main__":
    main()