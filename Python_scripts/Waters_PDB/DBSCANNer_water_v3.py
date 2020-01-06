import numpy as np
import matplotlib.pyplot as plt
#%matplotlib inline
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import DBSCAN
import os,re
from pymol import cmd
import matplotlib.cm as cm
import matplotlib.colors as colors

################### USER VARIABLES
my_eps = 0.9 # epsilon value for DBSCAN
my_n = 2 # n value for DBSCAN
occupancy_treshold = 3 # minimal occupancy for cluster to be considered
color_map = 'Reds'


############# Usage:
# run in pymol, while having structure loaded
# structure: all nuclear receptors aligned -> only best fitting structure kept -> waters superimposed on this one

file_array = []
working_dir = os.getcwd() + '/'
for r, d, f in os.walk(working_dir):
    for file in f:
        if '.pdb' in file:
            file_array.append(file)

if len(file_array) > 1:
    print('ERROR: found more than one file, will use: ' + file_array[0])
else:
    pass

water_array = []

current_pdb_open = open(file_array[0],'r')
current_pdb = current_pdb_open.read()
water_pattern = r'HETATM\s?\d+\s+\S+\s+HOH\s+\S+\s?\s?\S+\s+(\S+)\s+(\S+)\s+(\S+)'
watter_pattern_compile = re.compile(water_pattern)
water_matches = watter_pattern_compile.finditer(current_pdb)
for i1 in water_matches:
    sub_array = [float(i1.group(1)), float(i1.group(2)), float(i1.group(3))]
    water_array.append(sub_array)

water_pattern2 = r'HETATM\s?\d+\s+\S+\s+HOH\s+\S+\s\s\s\S+\s+(\S+)\s+(\S+)\s+(\S+)'
watter_pattern_compile2 = re.compile(water_pattern2)
water_matches2 = watter_pattern_compile2.finditer(current_pdb)
for i2 in water_matches2:
    sub_array = [float(i2.group(1)), float(i2.group(2)), float(i2.group(3))]
    water_array.append(sub_array)

water_pattern3 = r'HETATM\s?\d+\s+\S+\s+AHOH\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)'
watter_pattern_compile3 = re.compile(water_pattern3)
water_matches3 = watter_pattern_compile3.finditer(current_pdb)
for i3 in water_matches3:
    sub_array = [float(i3.group(1)), float(i3.group(2)), float(i3.group(3))]
    water_array.append(sub_array)

water_pattern4 = r'HETATM\s?\d+\s+\S+\s+BHOH\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)'
watter_pattern_compile4 = re.compile(water_pattern4)
water_matches4 = watter_pattern_compile4.finditer(current_pdb)
for i4 in water_matches4:
    sub_array = [float(i3.group(1)), float(i3.group(2)), float(i3.group(3))]
    water_array.append(sub_array)

print('Found ' + str(len(water_array)) + ' water molecules')

data = np.asarray(water_array)
fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(data[:,0], data[:,1], data[:,2], s=15)
ax.view_init(azim=200)
plt.show()

model = DBSCAN(eps=my_eps, min_samples=2)
model.fit_predict(data)
pred = model.fit_predict(data)

fig = plt.figure()
ax = Axes3D(fig)
ax.scatter(data[:,0], data[:,1], data[:,2], c=model.labels_, s=15)
ax.view_init(azim=200)
plt.show()

print("number of cluster found: {}".format(len(set(model.labels_))))
print('cluster for each point: ', model.labels_)
print('length?')
print(len(model.labels_))
print(len(set(model.labels_)))
cluster_outfile = open('clusters_out.tmp', 'w')

for i in range(0,((len(model.labels_)))):
    cur_coords = data[i]
    cur_clus = model.labels_[i]
    coords_flat = ''
    for aa in cur_coords:
        coords_flat = coords_flat + str(aa) + ' '

    string = str(coords_flat) + ' '  + str(cur_clus)
    print(string)
    cluster_outfile.write(string)
    cluster_outfile.write('\n')
cluster_outfile.close()

print('separating clusters...')
cluster_outfile_in_open = open('clusters_out.tmp', 'r')
cluster_outfile_in = cluster_outfile_in_open.read()

for cluster_nr in range(0,(len(set(model.labels_)))):
    cluster_individual_file_name = 'cluster_individual_' + str(cluster_nr) + '.ctmp'
    cluster_individual_file = open(cluster_individual_file_name,'w')

    c_pattern = r'(\S+\s+\S+\s+\S+)\s+(' + re.escape(str(cluster_nr)) + ')\s'
    c_pattern_compile = re.compile(c_pattern)
    c_matches = c_pattern_compile.finditer(cluster_outfile_in)
    for i5 in c_matches:
        cluster_individual_file.write(i5.group(1))
        cluster_individual_file.write('\n')

    cluster_individual_file.close()

ctmp_file_array = []
for r, d, f in os.walk(working_dir):
    for file in f:
        if '.ctmp' in file:
            ctmp_file_array.append(file)
removal_counter = 0
for ctmp_file in ctmp_file_array:
    ctmp_file_open = open(ctmp_file, 'r')
    ctmp_file_content = ctmp_file_open.read().split('\n')

    if len(ctmp_file_content) <= 1:
        print('removing file with content: ')
        print(ctmp_file_content)
        os.remove(ctmp_file)
        removal_counter = removal_counter + 1
    else:
        pass


print('calculating centroids and weigths')
cluster_middle_point_file = open('cluster_middle_points.res', 'w')

for cluster_file in range(0,(len(set(model.labels_))-removal_counter)):
    cluster_ind_filename = 'cluster_individual_' + str(cluster_file) + '.ctmp'
    current_cluster_file_open = open(cluster_ind_filename, 'r')
    current_cluster_file = current_cluster_file_open.read()

    cluster_x_coords = []
    cluster_y_coords = []
    cluster_z_coords = []

    coord_pattern = r'(\S+)\s+(\S+)\s+(\S+)'
    coord_pattern_compile = re.compile(coord_pattern)
    coord_matches = coord_pattern_compile.finditer(current_cluster_file)
    for i6 in coord_matches:

        cluster_x_coords.append(float(i6.group(1)))
        cluster_y_coords.append(float(i6.group(2)))
        cluster_z_coords.append(float(i6.group(3)))

    x_average = (sum(cluster_x_coords)) / (len(cluster_x_coords))
    y_average = (sum(cluster_y_coords)) / (len(cluster_y_coords))
    z_average = (sum(cluster_z_coords)) / (len(cluster_z_coords))
    print('centroid of cluster ' + str(cluster_file) + ' is located at ' + str(x_average) + ' ' + str(y_average) + ' ' + str(z_average) + ' and has a weight of ' + str(len(cluster_x_coords)))
    cluster_middle_point_file.write('centroid of cluster ' + str(cluster_file) + ' is located at ' + str(x_average) + ' ' + str(y_average) + ' ' + str(z_average) + ' and has a weight of ' + str(len(cluster_x_coords)))
    cluster_middle_point_file.write('\n')

cluster_middle_point_file.close()

remove_array1 = []
for r, d, f in os.walk(working_dir):
    for file in f:
        if '.ctmp' in file:
            remove_array1.append(file)
for remove_file in remove_array1:
    os.remove(remove_file)


middle_points_clean = open('middle_points_clean.res', 'w')

print('removing clusters with occupancy below 4...')

middle_points_dirty_open = open('cluster_middle_points.res','r')
middle_points_dirty = middle_points_dirty_open.read()

weights_array = [] 
middle_point_pattern = r'centroid\sof\scluster\s\d+\sis\slocated\sat\s\S+\s\S+\s\S+\sand\shas\sa\sweight\sof\s(\d+)'
middle_point_pattern_compile = re.compile(middle_point_pattern)
middle_point_matches = middle_point_pattern_compile.finditer(middle_points_dirty)
for i7 in middle_point_matches:
    print(i7.group(0)) 
    weights_array.append(int(i7.group(1)))
    print('the occupancy is: ' + str(i7.group(1))) 
    if int(i7.group(1)) >= occupancy_treshold:
        middle_points_clean.write(i7.group(0))
        middle_points_clean.write('\n')
middle_points_dirty_open.close()
middle_points_clean.close()

dummy_coords_file_open = open('middle_points_clean.res','r')
dummy_coords_lines = dummy_coords_file_open.read().split('\n')
dummy_coords_lines = dummy_coords_lines[:(len(dummy_coords_lines)-1)]

min_occupancy = occupancy_treshold
print('WEIGHTS ARRAY:')
print(weights_array)
max_occupancy = max(weights_array)

print('lower color threshold: ' + str(min_occupancy))
print('upper color treshold: ' + str(max_occupancy))

norm = colors.Normalize(vmin=min_occupancy, vmax=max_occupancy)
f2rgb = cm.ScalarMappable(norm=norm, cmap=cm.get_cmap(color_map))
def f2hex(f2rgb, f):
    rgb = (f2rgb.to_rgba(f)[:3])
    return rgb

pseudo_counter = 1
for coord_line in dummy_coords_lines:
    ps_coord_pattern = r'at\s(\S+)\s(\S+)\s(\S+)\sand\shas\sa\sweight\sof\s(\d+)'
    ps_coord_pattern_compile = re.compile(ps_coord_pattern)
    ps_coord_matches = ps_coord_pattern_compile.finditer(coord_line)
    for i8 in ps_coord_matches:

        pseudo_x = float(i8.group(1))
        pseudo_y = float(i8.group(2))
        pseudo_z = float(i8.group(3))
        weight = int(i8.group(4))

        v1 = f2hex(f2rgb, weight)
        v1_r = v1[0] * 255
        v1_g = v1[1] * 255
        v1_b = v1[2] * 255
        color_tupple = (v1_r, v1_g, v1_b)
        print('color tupple: ')
        print(v1_r, v1_g, v1_b)
        colname = 'col' + str(pseudo_counter)
        current_color = cmd.set_color(colname, color_tupple)

    pseudoatom_name = 'ps' + str(pseudo_counter)
    cmd.pseudoatom(pseudoatom_name, pos=[pseudo_x,pseudo_y,pseudo_z], color=colname)
    cmd.show('dots', pseudoatom_name)
    pseudo_counter = pseudo_counter + 1

print(str(len(dummy_coords_lines)))


