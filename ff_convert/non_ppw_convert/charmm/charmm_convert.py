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


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <pdbfile>")
        sys.exit(1)
    infile = sys.argv[1]
    step1_pdb = f"converted_{infile}"

    process_pdb(infile, step1_pdb)


