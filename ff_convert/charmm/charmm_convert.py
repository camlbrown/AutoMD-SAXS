import sys

def rename_glu(residue_lines):
    """Convert GLU to GLUP if HE2 present"""
    hydrogens = {line[12:16].strip() for line in residue_lines}
    return 'GLUP' if 'HE2' in hydrogens else 'GLU'


def rename_asp(residue_lines):
    """Convert ASP to ASPP if HD2 present"""
    hydrogens = {line[12:16].strip() for line in residue_lines}
    return 'ASPP' if 'HD2' in hydrogens else 'ASP'


def rename_his(residue_lines):
    """Determine protonation state: HSP (HIP), HSE, or HSD"""
    hydrogens = {line.split()[2] for line in residue_lines}
    if {'HE1','HE2','HD1','HD2'}.issubset(hydrogens):
        return 'HSP'
    if {'HE1','HE2','HD2'}.issubset(hydrogens):
        return 'HSE'
    if {'HD1','HD2','HE1'}.issubset(hydrogens):
        return 'HSD'
    return None


def process_pdb(infile, outfile):
    """
    Reads PDB, renames GLU->GLUP, ASP->ASPP, HIS->HSP/HSE/HSD, writes step1 PDB.
    """
    with open(infile) as f:
        lines = f.readlines()

    output = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('ATOM') and line[17:20].strip() in ('GLU','ASP','HIS'):
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
            elif res_name == 'ASP':
                new_name = rename_asp(residue_lines)
            else:
                new_name = rename_his(residue_lines)
            for l in residue_lines:
                if new_name:
                    output.append(f"{l[:17]}{new_name:<4}{l[21:]}")
                else:
                    output.append(l)
        else:
            output.append(line)
            i += 1

    with open(outfile, 'w') as f:
        f.writelines(output)


def convert_for_charmm(input_pdb, output_pdb):
    """
    Applies CHARMM-compatible atom-name mappings per residue_dict
    """
    residue_dict = {
        "ACE": {"1H": "HH31", "2H": "HH32", "3H": "HH33"},
        "GLU": {"HB2": "HB1 ", "HB3": "HB2 ", "HG2": "HG1 ", "HG3": "HG2 ", "HA2": " HA ", " H ": " HN", "H2": "HN  "},
        "VAL": {"HB2": " HB ", "HB3": " HB ", "HA2": " HA ", " H ": " HN"},
        "GLN": {"H1": " H  ", "H2": " H  ", "HB2": "HB1 ", "HB3": "HB2 ", "HG2": "HG1 ", "HG3": "HG2 ", "HA2": " HA ", " H ": " HN"},
        "LEU": {"HB2": "HB1 ", "HB3": "HB2 ", "H1": " H  ", "HXT": " H  ", "HA2": " HA ", " H ": " HN"},
        "SER": {"HB2": "HB1 ", "HB3": "HB2 ", "HA2": " HA ", " H ": " HN", "HG": "HG1 "},
        "ARG": {"HB2": "HB1 ", "HB3": "HB2 ", "HG2": "HG1 ", "HG3": "HG2 ", "HD2": "HD1 ", "HD3": "HD2 ", "HA2": " HA ", " H ": " HN"},
        "CYS": {"HB2": "HB1 ", "HB3": "HB2 ", "HA2": " HA ", " H ": " HN", "HG": "HG1 "},
        "TYR": {"HB2": "HB1 ", "HB3": "HB2 ", "HA2": " HA ", " H ": " HN"},
        "ILE": {"HG12": "HG11", "HG13": "HG12", "HD11": "HD1 ", "HD12": "HD2 ", "HD13": "HD3 ", "HA2": " HA ", "HB2": " HB ", " H ": " HN"},
        "TRP": {"HB2": "HB1 ", "HB3": "HB2 ", "HA2": " HA ", " H ": " HN"},
        "PRO": {"HD2": "HD1 ", "HD3": "HD2 ", "HG2": "HG1 ", "HG3": "HG2 ", "HB2": "HB1 ", "HB3": "HB2 ", "HA2": " HA ", " H ": " HN"},
        "LYS": {"HB2": "HB1 ", "HB3": "HB2 ", "HG2": "HG1 ", "HG3": "HG2 ", "HD2": "HD1 ", "HD3": "HD2 ", "HE2": "HE1 ", "HE3": "HE2 ", "HA2": " HA ", "H1": " H  ", " H ": " HN"},
        "ASN": {"HB2": "HB1 ", "HB3": "HB2 ", "H1": " H  ", "HD1": "HD21", "HD2": "HD22", "HA2": " HA ", " H ": " HN"},
        "ASP": {"HB2": "HB1 ", "HB3": "HB2 ", "H1": " H  ", "HA2": " HA ", " H ": " HN"},
        "MET": {"HB2": "HB1 ", "HB3": "HB2 ", "HG2": "HG1 ", "HG3": "HG2 ", "HA2": " HA ", " H ": " HN"},
        "NME": {"CA": "CH3 ", "1HA": "HH31", "2HA": "HH32", "3HA": "HH33", " H ": " HN"},
        "GLY": {"HA2": "HA1 ", "HA3": "HA2 ", " H ": " HN"},
        "PHE": {"H1": " H  ", "HB2": "HB1 ", "HB3": "HB2 ", "HA2": " HA ", " H ": " HN"},
        "THR": {"H1": " H  ", "H2": " H  ", "HA2": " HA ", " H ": " HN"},
        "ASPP": {"HB2": "HB1 ", "HB3": "HB2 ", "HA2": " HA ", "H2": " HN ", " H ": " HN"},
        "HIS": {"HB2": "HB1 ", "HB3": "HB2 ", "HA2": " HA ", " H ": " HN"},
        "HSE": {"HB2": "HB1 ", "HB3": "HB2 ", "HA2": " HA ", " H ": " HN"},
        "HSD": {"HB2": "HB1 ", "HB3": "HB2 ", "HA2": " HA ", " H ": " HN"},
        "GLUP": {"HB2": "HB1 ", "HB3": "HB2 ", "HG2": "HG1 ", "HG3": "HG2 ", "HA2": " HA ", "H2": " HN ", " H ": " HN"},
        "HSP": {"HB2": "HB1 ", "HB3": "HB2 ", "HA2": " HA ", " H ": " HN"},
        "ALA": {"H2": " H  ", "HA2": " HA ", " H ": " HN"}
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

    convert_for_charmm(step1_pdb, converted_pdb)

