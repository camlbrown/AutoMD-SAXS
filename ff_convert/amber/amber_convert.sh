#!/bin/bash
# Convert pdb to AMBER format

input_pdb=$1

sed -i '/^ANISOU/d' $input_pdb
sed -i '/NAG/d' $input_pdb
sed -i '/GOL/d' $input_pdb
sed -i '/DOD/d' $input_pdb
sed -i '/PEG/d' $input_pdb
sed -i '/PEG/d' $input_pdb
sed -i '/BAL/d' $input_pdb
sed -i '/EDO/d' $input_pdb
sed -i '/HOH/d' $input_pdb #future option to keep xtal waters
sed -i 's/NMA/NME/g' $input_pdb

python amber_convert.py $input_pdb
rm $input_pdb
rm step*
mv converted* GMX.pdb
