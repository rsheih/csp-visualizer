# csp-visualizer
A Python package for mapping chemical shift perturbations onto PDB/CIF structures for visualization in PyMOL and Chimera. 

### Notes 
- This tool replaces the B-factor column in the structure file with CSP values aligned by residue number
- The output structure can be directly visualized in PyMOL or Chimera using standard coloring by B-factor  
- - Do NOT use `~` (tilde) in file paths. Always use the full absolute path (e.g., `/Users/yourname/data/file.pdb`) instead of shortcuts like `~/data/file.pdb`.

### Prerequisite

- You must have **Python installed** on your system  
- If you do not have Python, download and install it here: https://www.python.org/downloads/

### How to Run the Program

1. **Download the files**  
   Download the entire `csp-visualizer-run` folder from this repository.

2. **Keep files together**  
   Make sure all files inside `csp-visualizer-run` stay in the same directory on your computer. Do not move individual files out of this folder.

3. **Open a terminal and navigate to the folder**  
   In your terminal, move into the folder using:
   ```bash
   cd /path/to/csp-visualizer-run

4. **Run the following command in the terminal:**
   ```bash
   bash run.sh 

5. **Provide the input files**  
   The program will prompt you to enter your input files. When prompted, use **absolute file paths** (the full path to the file on your system).
    Example:
   ```bash
   /Users/yourname/Documents/data/protein.pdb

Happy CSP visualizing! 🧬


