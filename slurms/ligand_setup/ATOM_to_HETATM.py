standard_residues = [
    "THR", "ACE", "GLU", "VAL", "GLN", "LEU", "SOL", "SER", "ARG", "HOH", "CYS", "TYR", "ILE", "TRP", "PRO", "LYS",
    "ASN", "ASP", "MET", "NME", "ACE", "GLY", "PHE", "HIS", "HIE", "HID", "HIP", "PHE", "ALA", "NHE", "CTHR", "CACE",
    "CGLU", "CVAL", "CGLN", "CLEU", "CSER", "CARG", "CCYS", "CCYX", "CYX", "NCYX", "CTYR", "CILE", "CTRP", "CPRO", "CLYS",
    "CASN", "CASP", "CMET", "CGLY", "CPHE", "CHIS", "CHIE", "HcID", "CHIP", "CPHE", "CALA", "GLH", "CYM", "CCYM", "NCYM",
    "HYP", "CHYP", "NHYP", "LYN", "NLYN", "ORN", "DAB", "ASH", "CGLY", "NGLH", "CASH", "NASH", "NTHR", "NGLU", "NVAL",
    "NGLN", "NLEU", "NSER", "NARG", "NCYS", "NTYR", "NILE", "NH2", "URE", "HO4", "NTRP", "NPRO", "NLYS", "NASN", "NASP",
    "NMET", "NGLY", "NPHE", "NHIS", "NHIE", "NHID", "NHIP", "NPHE", "NALA"
]

input_pdb_file = 'Complex_tidy.pdb'
output_pdb_file = 'GMX.pdb'

with open(input_pdb_file, 'r') as file:
    pdb_lines = file.readlines()

for i, line in enumerate(pdb_lines):
    if line.startswith('ATOM'):
        residue = line[17:20].strip()

        if residue not in standard_residues:
            pdb_lines[i] = line.replace('ATOM  ', 'HETATM', 1)

with open(output_pdb_file, 'w') as output_file:
    output_file.writelines(pdb_lines)

print(f"Updated PDB written to {output_pdb_file}")

