import os

ligand_dir = "ligands"
os.makedirs(ligand_dir, exist_ok=True)

with open("complex.pdb", "r") as pdb_file:
    protein_pdb = []
    ligand_pdb_dict = {}

    for line in pdb_file:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            residue_name = line[17:20].strip()
            residue_number = int(line[22:26])

            if residue_name not in ["THR", "ACE", "GLU", "VAL", "GLN", "LEU", "SOL", "SER", "ARG", "HOH", "CYS", "TYR", "ILE", "TRP", "PRO", "LYS", "ASN", "ASP", "MET", "NME", "ACE", "GLY", "PHE", "HIS", "HIE", "HID", "HIP", "PHE", "ALA", "NHE",  "CTHR", "CACE", "CGLU", "CVAL", "CGLN", "CLEU", "CSER", "CARG", "CCYS", "CCYX", "CYX", "NCYX", "CTYR", "CILE", "CTRP", "CPRO", "CLYS", "CASN", "CASP", "CMET", "CGLY", "CPHE", "CHIS", "CHIE", "HcID", "CHIP", "CPHE", "CALA", "GLH", "CYM", "CCYM", "NCYM", "HYP", "CHYP", "NHYP", "LYN", "NLYN", "ORN", "DAB",   "ASH", "CGLY", "NGLH", "CASH", "NASH", "NTHR", "NGLU", "NVAL", "NGLN", "NLEU", "NSER", "NARG", "NCYS", "NTYR", "NILE", "NH2", "URE", "HO4", "NTRP", "NPRO", "NLYS", "NASN", "NASP", "NMET", "NGLY", "NPHE", "NHIS", "NHIE", "NHID", "NHIP", "NPHE", "NALA", "MSE"]:
                # Check if the ligand with the same name and number already exists
                ligand_key = (residue_name, residue_number)
                if ligand_key not in ligand_pdb_dict:
                    ligand_pdb_dict[ligand_key] = []
                ligand_pdb_dict[ligand_key].append(line)
            else:
                protein_pdb.append(line)
        else:
            # Keep lines that are not ATOM or HETATM records (e.g., HEADER, TITLE)
            protein_pdb.append(line)

with open("protein.pdb", "w") as protein_file:
    protein_file.write("".join(protein_pdb))

for ligand_key, ligand_data in ligand_pdb_dict.items():
    ligand_residue_name, ligand_residue_number = ligand_key
    ligand_file_path = os.path.join(ligand_dir, f"ligand_{ligand_residue_number}.pdb")
    with open(ligand_file_path, "w") as ligand_file:
        ligand_file.write("".join(ligand_data))

