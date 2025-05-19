import re
import os
import subprocess
import sys

# Renaming functions for different residues based on hydrogen content
def rename_glu(residue_lines):
    hydrogens = {line.split()[2] for line in residue_lines}
    return 'GLH' if 'HE2' in hydrogens else 'GLU'


def rename_his(residue_lines):
    hydrogens = {line.split()[2] for line in residue_lines}
    if {'HE1','HE2','HD1','HD2'}.issubset(hydrogens):
        return 'HIP'
    if {'HE1','HE2','HD2'}.issubset(hydrogens):
        return 'HIE'
    if {'HD1','HD2','HE1'}.issubset(hydrogens):
        return 'HID'
    return None


def rename_asp(residue_lines):
    hydrogens = {line.split()[2] for line in residue_lines}
    if 'H2' in hydrogens:
        return 'NASP'
    if 'HD2' in hydrogens:
        return 'ASH'
    return 'ASP'


def remove_whitespace(line, start=20, end=27):
    # Remove the first whitespace in the specified region
    return line[:start] + line[start:end].replace(' ', '', 1) + line[end:]


def process_pdb(infile, outfile):
    """
    Reads a PDB file, renames GLU, HIS, and ASP residues based on hydrogen content,
    and writes the result to outfile.
    """
    with open(infile) as f:
        lines = f.readlines()

    output = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('ATOM') and line[17:20] in ('GLU','HIS','ASP'):
            res_name = line[17:20]
            res_id = line[22:26].strip()
            residue_lines = []
            while i < len(lines) and lines[i].startswith('ATOM') \
                  and lines[i][17:20] == res_name \
                  and lines[i][22:26].strip() == res_id:
                residue_lines.append(lines[i])
                i += 1
            if res_name == 'GLU':
                new_name = rename_glu(residue_lines)
            elif res_name == 'HIS':
                new_name = rename_his(residue_lines)
            else:
                new_name = rename_asp(residue_lines)
            for l in residue_lines:
                if new_name:
                    if res_name == 'ASP' and len(new_name) != 3:
                        l2 = remove_whitespace(l)
                        output.append(l2[:17] + new_name + l2[20:])
                    else:
                        output.append(l[:17] + new_name + l[20:])
                else:
                    output.append(l)
        else:
            output.append(line)
            i += 1

    with open(outfile, 'w') as f:
        f.writelines(output)


def convert_for_gmx(input_pdb, output_pdb):
    """
    Applies atom-name mappings for Gromacs compatibility based on residue_dict,
    writing converted PDB to output_pdb.
    """
    residue_dict = {
        "ACE": {"1H": "HH31", "2H": "HH32", "3H": "HH33"},
        "GLU": {"HB2": "HB1 ", "HB3": "HB2 ", "HG2": "HG1 ", "HG3": "HG2 "},
        "GLU": {"HB2": "HB1 ", "HB3": "HB2 ", "HG2": "HG1 ", "HG3": "HG2 "},
        "VAL": {"HB2": " HB ", "HB3": " HB "},
        "GLN": {"HB2": "HB1 ", "HB3": "HB2 ", "HG2": "HG1 ", "HG3": "HG2 "},
        "LEU": {"HB2": "HB1 ", "HB3": "HB2 "},
        "SER": {"HB2": "HB1 ", "HB3": "HB2 "},
        "ARG": {"HB2": "HB1 ", "HB3": "HB2 ", "HG2": "HG1 ", "HG3": "HG2 ", "HD2": "HD1 ", "HD3": "HD2 "},
        "CYS": {"HB2": "HB1 ", "HB3": "HB2 "},
        "TYR": {"HB2": "HB1 ", "HB3": "HB2 "},
        "ILE": {"HG12": "HG11", "HG13": "HG12", "HD11": "HD1 ", "HD12": "HD2 ", "HD13": "HD3 "},
        "TRP": {"HB2": "HB1 ", "HB3": "HB2 "},
        "PRO": {"HD2": "HD1 ", "HD3": "HD2 ", "HG2": "HG1 ", "HG3": "HG2 ", "HB2": "HB1 ", "HB3": "HB2 ", "H2": "H1", "H3": "H2"},
        "LYS": {"HB2": "HB1 ", "HB3": "HB2 ", "HG2": "HG1 ", "HG3": "HG2 ", "HD2": "HD1 ", "HD3": "HD2 ", "HE2": "HE1 ", "HE3": "HE2 "},
        "ASN": {"HB2": "HB1 ", "HB3": "HB2 ", "H1": " H  ", "HD1": "HD21", "HD2": "HD22", "HA2": " HA "},
        "ASP": {"HB2": "HB1 ", "HB3": "HB2 "},
        "MET": {"HB2": "HB1 ", "HB3": "HB2 ", "HG2": "HG1 ", "HG3": "HG2 "},
        "NME": {"CA": "CH3 ", "1HA": "HH31", "2HA": "HH32", "3HA": "HH33"},
        "NMA": {"CA": "CH3 ", "1HA": "HH31", "2HA": "HH32", "3HA": "HH33"},
        "GLY": {"HA2": "HA1 ", "HA3": "HA2 "},
        "PHE": {"HB2": "HB1 ", "HB3": "HB2 "},
#        "THR": {},
        "ASH": {"HB2": "HB1 ", "HB3": "HB2 "},
        "HIS": {"HB2": "HB1 ", "HB3": "HB2 "},
        "HIE": {"HB2": "HB1 ", "HB3": "HB2 "},
        "HID": {"HB2": "HB1 ", "HB3": "HB2 "},
        "GLH": {"HB2": "HB1 ", "HB3": "HB2 ", "HG2": "HG1 ", "HG3": "HG2 "},
        "HIP": {"HB2": "HB1 ", "HB3": "HB2 "},
        "NASP": {"HB2": "HB1 ", "HB3": "HB2 "},
#        "ALA": {}
    }


    with open(input_pdb) as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            if line.startswith(('ATOM','HETATM')):
                resname = line[17:21].strip()
                atom_name = line[12:16].strip()
                if resname in residue_dict and atom_name in residue_dict[resname]:
                    new_atom = residue_dict[resname][atom_name].ljust(4)
                    line = line[:12] + new_atom + line[16:]
            f_out.write(line)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <pdbfile>")
        sys.exit(1)
    infile = sys.argv[1]
    step1_pdb = f"step1_{infile}"
    converted_pdb = f"converted_{infile}"

    process_pdb(infile, step1_pdb)

    convert_for_gmx(step1_pdb, converted_pdb)

