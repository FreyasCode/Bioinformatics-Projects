import random

amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
for i in range(20):
    peptide_chain = ''.join(random.choice(amino_acids) for _ in range(random.randint(20, 40)))
    print(peptide_chain)

# trimer_dict = {}
# trimer_count = len(peptide_chain) // 3
# length = len(peptide_chain)
# for i in range(0, len(peptide_chain) - 2, 3): # TODO: step every 1 or 3? 
#     # Extract a trimer
#     if i+3 > length:
#         break
#     trimer = peptide_chain[i:i+3]
#     # Check if the trimer is already in the dictionary
#     if trimer in trimer_dict:
#         # Increment the count if it is
#         trimer_dict[trimer] += 1
#     else:
#         # Add the trimer to the dictionary with a count of 1 if it's not
#         trimer_dict[trimer] = 1

# p = {key: value / trimer_count for key, value in trimer_dict.items()}
# print(p)