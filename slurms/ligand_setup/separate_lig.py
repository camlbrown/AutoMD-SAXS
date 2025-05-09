import numpy

with open("complex.pdb", "r") as pdb_file:
    protein_pdb = []
    ligand_pdb = []
    ligand_section = False

    for line in pdb_file:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            residue_name = line[17:20].strip()
            
            if residue_name.isalpha() and len(residue_name) == 3:
                # Check if the residue name is not a standard amino acid (e.g., ALA, GLY)
                if residue_name not in ["THR", "ACE", "GLU", "VAL", "GLN", "LEU", "SOL", "SER", "ARG", "HOH", "CYS", "TYR", "ILE", "TRP", "PRO", "LYS", "ASN", "ASP", "MET", "NME", "ACE", "GLY", "PHE", "HIS", "HIE", "HID", "HIP", "PHE", "ALA", "NHE",  "CTHR", "CACE", "CGLU", "CVAL", "CGLN", "CLEU", "CSER", "CARG", "CCYS", "CCYX", "CYX", "NCYX", "CTYR", "CILE", "CTRP", "CPRO", "CLYS", "CASN", "CASP", "CMET", "CGLY", "CPHE", "CHIS", "CHIE", "HcID", "CHIP", "CPHE", "CALA", "GLH", "CYM", "CCYM", "NCYM", "HYP", "CHYP", "NHYP", "LYN", "NLYN", "ORN", "DAB",   "ASH", "CGLY", "NGLH", "CASH", "NASH", "NTHR", "NGLU", "NVAL", "NGLN", "NLEU", "NSER", "NARG", "NCYS", "NTYR", "NILE", "NH2", "URE", "HO4", "NTRP", "NPRO", "NLYS", "NASN", "NASP", "NMET", "NGLY", "NPHE", "NHIS", "NHIE", "NHID", "NHIP", "NPHE", "NALA"]:  
                    ligand_section = True
                    ligand_pdb.append(line)
                else:
                    protein_pdb.append(line)
            else:
                # Two-letter code may indicate an ion (e.g., PB)
                if len(residue_name) == 2:
                    protein_pdb.append(line)
                else:
                    ligand_section= True
                    ligand_pdb.append(line)
        else:
            # Keep lines that are not ATOM or HETATM records (e.g., HEADER, TITLE)
            protein_pdb.append(line)

# Write protein and ligand PDB files
with open("protein.pdb", "w") as protein_file:
    protein_file.write("".join(protein_pdb))

with open("ligand.pdb", "w") as ligand_file:
    ligand_file.write("".join(ligand_pdb))

