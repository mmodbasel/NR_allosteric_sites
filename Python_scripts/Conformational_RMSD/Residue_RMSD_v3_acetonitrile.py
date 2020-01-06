import numpy as np
import re
import pandas as pd
import os

"""
Calculates RMSD for a selection of residues between two input files

EXAMPLE INPUT:
ATOM      1  N   CYS A 669       9.744  10.513 -18.120  1.00 71.41           N1+

OUTPUT:
Q727: 2.1 A
A728  2.3 A
E729  3.2 A
...

IMPORTANT:
-expects residues to have same order of atoms and only account for atom types C, N, O , S
-needs list file to know which residues to be processed
"""

### USER DEFINED VARIABLES

reference_pdb_filename = 'AR.pdb'
analyte_pdb_filename = 'AR_acetonitrile.pdb'
solvent = 'acetonitrile'
receptor = 'AR'
site = 'AF2'
chain_ident = 'A'


### EXECUTION

output_filename = 'RMSD_' + receptor + '_' + solvent + '_res.out'
output_file = open(output_filename,'w')
index_dict = {'AR':1, 'ERa':2, 'ERB':3, 'GR':4, 'MR':5, 'PR':6, 'TRa':7, 'TRB':8}
resnr_array = []
list_name = receptor + '_' + site + '.lst'
list_open = open(list_name,'r')
list_lines = list_open.read().split('\n')
list_lines = list_lines[:-1]
for list_line in list_lines:
        number = (str(list_line))[1:]
        resnr_array.append(int(number))
AR_resnr_array = []
AR_list_file = receptor + '_' + site + '.lst'
AR_list_open = open(AR_list_file,'r')
AR_list_lines = AR_list_open.read().split('\n')
AR_list_lines = AR_list_lines[:-1]
for AR_list_line in AR_list_lines:
        AR_number = AR_list_line[1:]
        AR_resnr_array.append(int(AR_number))
csv_name = 'data_empty_' + site + '_' + solvent + '.csv'
data_file = pd.read_csv(csv_name)
for resnr in resnr_array:
        if resnr == 'nan':
                pass
        reference_file_open = open(reference_pdb_filename,'r')
        reference_file = reference_file_open.read()
        analyte_file_open = open(analyte_pdb_filename, 'r')
        analyte_file = analyte_file_open.read()
        coords_reference = []
        coords_analyte = []
        residue_pattern = r'ATOM\s+\S+\s+[CNOS]\S+\s+\S+\s+' + re.escape(chain_ident) + r'\s+' + re.escape(str(resnr)) + r'\s+(\S+)\s+(\S+)\s+(\S+)'
        residue_pattern_c = re.compile(residue_pattern)
        residue_matches_ana = residue_pattern_c.finditer(analyte_file)
        residue_matches_ref = residue_pattern_c.finditer(reference_file)
        for ref_match in residue_matches_ref:
                print(ref_match.group(0))
                res_coords = [ref_match.group(1), ref_match.group(2), ref_match.group(3)]
                coords_reference.append(res_coords)
        for ana_match in residue_matches_ana:
                print(ana_match.group(0))
                res_coords = [ana_match.group(1), ana_match.group(2), ana_match.group(3)]
                coords_analyte.append(res_coords)
        if len(coords_reference) != len(coords_analyte):
                print('FATAL ERROR: unequal amount of heavy atoms read from files')
        else:
                pass
        distances_atoms_current_residue = []
        for iterator in range(0,len(coords_reference)):
                ref_coord_array = coords_reference[iterator]
                ana_coord_array = coords_analyte[iterator]
                distance_vector_x = float(ref_coord_array[0]) - float(ana_coord_array[0])
                distance_vector_y = float(ref_coord_array[1]) - float(ana_coord_array[1])
                distance_vector_z = float(ref_coord_array[2]) - float(ana_coord_array[2])
                distance = np.sqrt( (distance_vector_x*distance_vector_x) + (distance_vector_y*distance_vector_y) + (distance_vector_z*distance_vector_z) )
                distances_atoms_current_residue.append(float(distance))
        residue_rmsd = sum(distances_atoms_current_residue) / len(distances_atoms_current_residue)
        print('RMSD is ' + str(residue_rmsd))
        output_file.write(str(resnr) + ': ' + str(residue_rmsd))
        output_file.write('\n')
        if receptor != 'AR':   # if receptor is not AR, residue numbers need to be translated back to the ones from AR
                current_index = (resnr_array.index(resnr))
                AR_resnr = AR_resnr_array[current_index]
                print('translating residue number back to AR numbering: ' + str(resnr) + '--> ' + str(AR_resnr))
        else:
                AR_resnr = resnr
        receptor_index = index_dict[receptor] -1
        print('receptor index: ' + str(receptor_index))
        data_file.at[receptor_index, str(AR_resnr)] = residue_rmsd
os.remove(csv_name)
data_file.to_csv(csv_name, index=False)


