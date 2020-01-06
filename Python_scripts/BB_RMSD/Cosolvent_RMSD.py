import os
import numpy as np
import subprocess
import shutil
import stat

# user variables

atoms_of_interest = ['N', 'CA', 'C', 'O'] #backbone array
def pdbReader(pdbDir):
    # sequentially reads PDB files in given directory and puts them into file-name-array
    pdb_filearray = []
    for r, d, f in os.walk(pdbDir):
        for file in f:
            if '.pdb' in file:
                pdb_filearray.append(file)
    pdb_filearray.sort()
    return pdb_filearray

def vectorCalc(x1,y1,z1,x2,y2,z2):
    # calculates vector length between two points (x1, y1, z1) and (x2, y2, z2)
    dist_vector_x = float(x2) - float(x1)
    dist_vector_y = float(y2) - float(y1)
    dist_vector_z = float(z2) - float(z1)

    dist_vector_sq_x = dist_vector_x*dist_vector_x
    dist_vector_sq_y = dist_vector_y*dist_vector_y
    dist_vector_sq_z = dist_vector_z*dist_vector_z

    dist = np.sqrt(dist_vector_sq_x + dist_vector_sq_y + dist_vector_sq_z)

    print('distance: ' + str(dist))
    return dist

if __name__ == '__main__':
    # reading all PDB Files in dir and subdir
    working_directory = os.getcwd() + '/'
    pdb_file_list = pdbReader(working_directory)

    #opening output file for RMSD result
    result_rmsd_file = open('result_RMSD.txt', 'w')
    result_rmsd_file.write('0')
    result_rmsd_file.write('\n')

    # reading first reference frame
    reference_filename = pdb_file_list[0]
    reference_open = open(reference_filename,'r')
    reference_file_content = reference_open.read().split('\n')
    reference_file_content = reference_file_content[:-1]    # removing last (emtpy) line
    if 'CRYST' in reference_file_content[0]:    # cutting crystal line if present
        reference_file_content[1:]
    print('Read the following PDB Files: ')
    print(pdb_file_list)
    
    # loop over remaining files
    for pdb_filename in pdb_file_list:
        current_frame_distances = []
        if pdb_filename != pdb_file_list[0]:

            # rmsd comparison

            rmsd_file_open = open(pdb_filename)
            rmsd_file_c = rmsd_file_open.read().split('\n')
            rmsd_file_c = rmsd_file_c[:-1]     # removing last (emtpy) line
            if 'CRYST' in rmsd_file_c[0]:
                rmsd_file_c[1:]
            if len(rmsd_file_c) != len(reference_file_content):
                print('ERROR: files do not have same length')
                break
            for line in rmsd_file_c:
                resname = line[17:20]
                #print('resname: ' + resname)
                if not resname == 'WAT':
                    if not 'ANIS' in line:
                        if not 'HEADER' in line:
                            if not 'TITLE' in line:
                                if not 'REMARK' in line:
                                    if not 'HETATM' in line:
                                        atom_name = line[13:17].strip()
                                        if atom_name in atoms_of_interest:
                                            current_line_index = rmsd_file_c.index(line) # line index

                                            #DEBUG & CHECK
                                            print('LINE INDEX: ' + str(current_line_index))
                                            print('reference:')
                                            print(reference_file_content[current_line_index])
                                            print('compare line: ')
                                            print(line)
                                            ref_line = reference_file_content[current_line_index]
                                            ref_x = ref_line[30:38].strip()
                                            print('ref_x: ' + str(ref_x))
                                            ref_y = ref_line[38:46].strip()
                                            print('ref_y: ' + str(ref_y))
                                            ref_z = ref_line[46:54].strip()
                                            print('ref_z: ' + str(ref_z))
                                            comp_x = line[30:38].strip()
                                            print('comp_x: ' + str(comp_x))
                                            comp_y = line[38:46].strip()
                                            print('comp_y: ' + str(comp_y))
                                            comp_z = line[46:54].strip()
                                            print('comp_z: ' + str(comp_z))
                                            # calculating distance for current line (matching BB criteria)
                                            atom_dist = vectorCalc(float(ref_x), float(ref_y), float(ref_z), float(comp_x), float(comp_y), float(comp_z))
                                            current_frame_distances.append(float(atom_dist))
            current_frame_rmsd = sum(current_frame_distances) / len(current_frame_distances)
            result_rmsd_file.write(str(np.round(current_frame_rmsd,decimals=2)) + ' ' + pdb_filename)
            result_rmsd_file.write('\n')
