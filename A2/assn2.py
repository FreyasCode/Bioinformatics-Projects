import json
import sys
import math
from joblib import load

# loading p values from json file
filename = "probabilities.json"
with open(filename, 'r') as file:
    loaded_dict = json.load(file)
q = 0.000125 # same probability was obtained empirically
amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
amino_acid_dict = {
 'A': 0,
 'R': 1,
 'N': 2,
 'D': 3,
 'C': 4,
 'E': 5,
 'Q': 6,
 'G': 7,
 'H': 8,
 'I': 9,
 'L': 10,
 'K': 11,
 'M': 12,
 'F': 13,
 'P': 14,
 'S': 15,
 'T': 16,
 'W': 17,
 'Y': 18,
 'V': 19
}
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
def calculate_score2(peptide, model):
    encoding = [0 for i in range(20)]
    for p in peptide:
      encoding[amino_acid_dict[p]] += 1
    return round(model.predict_proba([encoding])[:,1][0],4)

# __main__
peptides = []
score1 = []
score2 = []
model = load('trained_model.joblib')

with open(sys.argv[2], 'r') as file:
    for line in file:
        peptides.append(line.rstrip('\n'))
        score1.append(calculate_score1(line.rstrip('\n'), loaded_dict))
        score2.append(calculate_score2(line.rstrip('\n'), model))

with open(sys.argv[4], "w") as file:
    for i in range(len(peptides)):
        file.write(str(score1[i]) + "\t" + str(score2[i]) + "\t" + peptides[i] + "\n")
    file.close()
