from Bio import SeqIO
import json

# Path to your FASTA file
fasta_file = "uniprot_sprot.fasta"

# Initialize an empty list to store sequences
sequences = []
trimer_count = 0

# Use SeqIO.parse to read the FASTA file
for record in SeqIO.parse(fasta_file, "fasta"):
    sequences.append(str(record.seq))
    trimer_count += len(record.seq) // 3

# Now, 'sequences' contains all the sequences from the FASTA file
# print(sequences[0:5])
# print("done!")

trimer_dict = {}
for s in sequences:
    length = len(s)
# Iterate through the string in steps of 3
    for i in range(0, len(s) - 2, 3): # TODO: step every 1 or 3? 
        # Extract a trimer
        if i+3 > length:
            break
        trimer = s[i:i+3]
        # Check if the trimer is already in the dictionary
        if trimer in trimer_dict:
            # Increment the count if it is
            trimer_dict[trimer] += 1
        else:
            # Add the trimer to the dictionary with a count of 1 if it's not
            trimer_dict[trimer] = 1

    # Print the dictionary to see the trimers and their frequencies
# print(trimer_count)

p = {key: value / trimer_count for key, value in trimer_dict.items()}
q = 1 / (20**3)

# Convert the dictionary to a JSON string
filename = 'probabilities.json'

# Writing the dictionary to a JSON file
with open(filename, 'w') as file:
    json.dump(p, file, indent=4)