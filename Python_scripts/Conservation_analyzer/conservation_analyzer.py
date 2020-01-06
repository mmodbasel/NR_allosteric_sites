import os
import pandas
import matplotlib.pyplot as plt
import seaborn as sns

""" 
Takes for a receptor site (either AF2 or BF3) that contains residue codes for each position of the site
L L L V M L I I (AR, ERa, ERB, ... TRB) for 1 residue position
V V V I I L T T
V I I V V V V V
K K K K K K K K
V L L L L I C C
Q Q Q Q Q Q Q Q
...

Fills table with conservation degree based on Clustal W classification
- needs csv file called "results_conservation.csv" 
"""
### USER VARIABLES
file_name = 'AF2_Residues.txt'
csv_file_name = 'results_conservation.csv'
rec_array = ['AR', 'ERa', 'ERB', 'GR', 'MR', 'PR', 'TRa', 'TRB']

def resCompare(residue1, residue2):
    # accepts two residues in single letter amino acid code as strings and compares them accoring to set criteria
    global res1_result
    blue_residues = ["A", "I", "L", "M", "F", "W", "V", "C"]
    red_residues = ["K", "R"]
    magenta_residues = ["E", "D"]
    green_residues = ["N", "Q", "S", "T"]
    proline = ["P"]
    aromatic = ['H', 'Y']
    glycine = ['G']
    if residue1.upper() in blue_residues:
        res1_result = 'blue'
    elif residue1.upper() in red_residues:
        res1_result = 'red'
    elif residue1.upper() in magenta_residues:
        res1_result = 'magenta'
    elif residue1.upper() in green_residues:
        res1_result = 'green'
    elif residue1.upper() in proline:
        res1_result = 'proline'
    elif residue1.upper() in aromatic:
        res1_result = 'aromatic'
    elif residue1.upper() in glycine:
        res1_result = 'glycine'
    if residue2.upper() in blue_residues:
        res2_result = 'blue'
    elif residue2.upper() in red_residues:
        res2_result = 'red'
    elif residue2.upper() in magenta_residues:
        res2_result = 'magenta'
    elif residue2.upper() in green_residues:
        res2_result = 'green'
    elif residue2.upper() in proline:
        res2_result = 'proline'
    elif residue2.upper() in aromatic:
        res2_result = 'aromatic'
    elif residue2.upper() in glycine:
        res2_result = 'glycine'
    else:
        return "FATAL ERROR: Compare function does not get a value"
    if residue1 == residue2:
        result_value = 1.0
    elif res1_result == res2_result and not residue1 == residue2:
        result_value = 0.5
    else:
        result_value = 0
    return result_value

def linePositionConverter(pos):
    rec_array = ['AR', 'ERa', 'ERB', 'GR', 'MR', 'PR', 'TRa', 'TRB']
    receptor_name = rec_array[pos]
    return receptor_name
working_dir = os.getcwd() + '/'
csv_file = pandas.read_csv('results_conservation.csv')
print('read csv file: ')
print(csv_file.head())
file_open = open(file_name,'r')
line_array = file_open.read().split('\n')
line_array = line_array[:-1]
number_of_residues = len(line_array)
sanity_1 = 0
for line in line_array:
    line = line.replace(' ', '')
    print(line)
    for rec_num in range(0,8):      # loop for receptors: 0 = AR, 1 = ERa, ...
        for line_position in range(0,8):
            rec_index = str(linePositionConverter(line_position))
            print('comparing: ' + line[rec_num] + ' and ' + line[line_position])
            print('would add the result at line ' + str(rec_num) + ' row ' + str(line_position) + ' into csv file')
            compare_result = resCompare(line[rec_num],line[line_position])
            print('result: ' + str(compare_result))
            old_value = float(csv_file.loc[int(rec_num), rec_index])
            new_value = old_value + float(compare_result)
            csv_file.at[int(rec_num), rec_index] = new_value
            if rec_num == 0 and line_position == 0:
                sanity_1 = sanity_1 + 1
            else:
                pass
if int(sanity_1) != int(number_of_residues):
    print('ERROR: number of residues to compare does not match number of times inner for loop executed')
    print(str(number_of_residues) + ' vs ' + str(sanity_1))
else:
    print('sanity check 1 passed!! ')
new_csv_path = working_dir + 'data_filled.csv'
csv_file.to_csv(new_csv_path)
zeilen = [0,1,2,3,4,5,6,7]
spalten = ['AR', 'ERa', 'ERB', 'GR', 'MR', 'PR', 'TRa', 'TRB']
for zeile in zeilen:
    for spalte in spalten:
        current_value = int(csv_file.loc[zeile,spalte])
        percentage = int(current_value) / int(number_of_residues)
        percentage = percentage * 100
        print('getting percentage for value ' + str(current_value))
        print('percentage is: ' + str(percentage))
        csv_file.at[zeile,spalte] = percentage
percentage_csv_path = working_dir + 'data_filled_percentage.csv'
csv_file.to_csv(percentage_csv_path)
sns.heatmap(csv_file, annot=True, cmap='Blues',fmt='g')
savepath = os.getcwd() + '/' + 'conserv1.png'
plt.savefig(savepath, dpi=1200)





