#!/usr/bin/env python3

# SCRIPT USED WHEN USER HAS PDB ID

#####################
# CONVERT XLSX FILE #
#####################
import os
import pandas as pd 
import time
print('imported')
from biopandas.pdb import PandasPdb
import sys
print(sys.executable)
# test: 2F1W

def get_structure(pdb_id, nmr_data): 
    'Load user inputted structure ID and NMR data file'
    df = pd.read_excel(nmr_data,  sheet_name=1) # make sure to ask user to give header
    print('loaded excel')

    nmr_data = df.rename(columns={"Unnamed: 0" : "Residue", "1to4" : "CSP"}) # rename column headers

    structure=PandasPdb().fetch_pdb(pdb_id) # fetch a copy of PDB file  
    if structure is None:
        print("No valid structure loaded.") 
    else:
        print('loaded structure')

    return structure, nmr_data

def substitute_b_factor(structure_id, nmr_data, save_dir=None):  
    '''Substitutes b-factor column with chemical shift perturbation data; will output PDB file'''
    structure, nmr_data = get_structure(structure_id, nmr_data) 
    # csp_df = pd.read_excel(csp_file) # excel file with chem shift info used to substitute b factor values

    if structure is None:
        print("No valid structure loaded. Exiting.")
        return None
    
    # rename column headers
    nmr_data.columns = ['Residue', 'CSP']
    atom_df=structure.df['ATOM'].copy() # get data with atom key
    # print(structure.df['ATOM'].columns) # check key names in file 

    # print(atom_df[['residue_number', 'residue_name', 'b_factor']].head(10)) # check before refactoring
    nmr_data = dict(zip(nmr_data['Residue'], nmr_data['CSP'])) # convert csp pandas df to file to prepare for merge
    atom_df['b_factor'] = atom_df['residue_number'].map(nmr_data).fillna(0.0) # map values to residues

    # check after merge
    # check_df = atom_df[['residue_number', 'residue_name', 'b_factor']] 
    # print(check_df.head(20))

    # put results back into data frame
    structure.df['ATOM'] = atom_df 

    # define save directory
    if save_dir is None:
        save_dir = os.getcwd()
    
    out_filename = f"demo_{structure_id}_b-factor.pdb"
    out_path = os.path.join(save_dir, out_filename)
    return structure, out_path
 
def return_file(structure, out_path):
    return structure.to_pdb(path=out_path, records=['ATOM'], gz=False), f'File is located at {out_path}'

def prompt_user(): 
    save_dir = None
    pdb_id = input("Input a PDB ID: ")
    nmr_data = input("Input the file path to your data file: ")
    save_dir = input("Input directory to save output (or press enter for current dir): ")
    save_dir = save_dir if save_dir.strip() else None
    return pdb_id, nmr_data, save_dir

def main(): 
    pdb_id, nmr_data, save_dir = prompt_user()
    structure, out_path = substitute_b_factor(pdb_id, nmr_data, save_dir) 
    new_file=return_file(structure, out_path)
    return new_file

if __name__ == "__main__":
    main()