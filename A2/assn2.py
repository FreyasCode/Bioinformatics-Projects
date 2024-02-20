import json
import sys
import math

# loading p values from json file
filename = "probabilities.json"
with open(filename, 'r') as file:
    loaded_dict = json.load(file)
q = 0.000125
amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# # Displaying the loaded dictionary
# print('Loaded dictionary from JSON file:')
# print(loaded_dict)

# Function to calculate score1
def calculate_score1(peptide, probabilities):
    score = 0
    for i in range(0, len(peptide)-2, 3): # TODO: check lengths not divisible by 3
        k_mer = peptide[i:i+3]
        if k_mer in probabilities:
            score += math.log2(probabilities[k_mer] / q)
        else:
            score += 0
    return round(score, 2)

# Function to calculate score2

# __main__
peptides = []
score1 = []
score2 = []
with open(sys.argv[2], 'r') as file:
    for line in file:
        peptides.append(line.rstrip('\n'))
        score1.append(calculate_score1(line.rstrip('\n'), loaded_dict))

with open(sys.argv[4], "w") as file:
    for i in range(len(peptides)):
        file.write(str(score1[i]) + "\t" + peptides[i] + "\n")
    file.close()

print(score1)

