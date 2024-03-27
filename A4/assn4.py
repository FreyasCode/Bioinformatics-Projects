from pyteomics import mgf, mass
import numpy as np
import sys
import math

# calculate spectrum mass (not m/z)
def calculate_spectrum_mass(pepmass, charge):
    return (pepmass - 1.0073) * charge

# convert mass to m/z for peptide in database
amino_acid_masses = {
    'A': 71.03711,  # Alanine
    'R': 156.10111, # Arginine
    'N': 114.04293, # Asparagine
    'D': 115.02694, # Aspartic Acid
    'C': 103.00919, # Cysteine
    'E': 129.04259, # Glutamic Acid
    'Q': 128.05858, # Glutamine
    'G': 57.02146,  # Glycine
    'H': 137.05891, # Histidine
    'I': 113.08406, # Isoleucine
    'L': 113.08406, # Leucine
    'K': 128.09496, # Lysine
    'M': 131.04049, # Methionine
    'F': 147.06841, # Phenylalanine
    'P': 97.05276,  # Proline
    'S': 87.03203,  # Serine
    'T': 101.04768, # Threonine
    'W': 186.07931, # Tryptophan
    'Y': 163.06333, # Tyrosine
    'V': 99.06841   # Valine
}

def read_peptides(filepath):
    peptides = {}
    with open(filepath, 'r') as file:
        for line in file:
            peptide, mass = line.strip().split('\t')
            peptides[peptide] = float(mass)
    return peptides

def generate_theoretical_spectrum(peptide):
    ions = []
    for i in range(1, len(peptide)):
        for charge in range(1, 3): 
            b_ion = mass.fast_mass(peptide[:i], ion_type='b', charge=charge)
            y_ion = mass.fast_mass(peptide[i:], ion_type='y', charge=charge) # TODO: check charge
            ions.extend([b_ion, y_ion])
    return np.array(ions)

def discretize_spectrum(mz_array, intensities, min_mz, max_mz, bin_width):
    num_bins = int(np.ceil((max_mz - min_mz) / bin_width))
    discrete_spectrum = np.zeros(num_bins)

    for i in range(len(mz_array)):
        bin_index = int((mz_array[i] - min_mz) / bin_width)
        if intensities is not None: 
            discrete_spectrum[bin_index] += math.sqrt(intensities[i]) # TODO: sqrt(intensities[i])
        else:
            discrete_spectrum[bin_index] += 1

    return discrete_spectrum

def cosine_similarity(v1, v2):
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    return np.dot(v1 / norm_v1, v2 / norm_v2) if norm_v1 and norm_v2 else 0

# check if mass difference is within tolerance
# if yes, evaluate score
def score_peptide_spectrum_match_intensity(peptide, spectrum, min_mz, max_mz, bin_size=0.5): # TODO: best value
    undiscretized_theoretical_spectrum = generate_theoretical_spectrum(peptide)
    theoretical_spectrum_vector = discretize_spectrum(undiscretized_theoretical_spectrum, None, min_mz, max_mz, bin_size)
    
    actual_spectrum_mz = np.array(spectrum['m/z array'])
    actual_spectrum_intensities = np.array(spectrum['intensity array'])
    actual_spectrum_vector = discretize_spectrum(actual_spectrum_mz, actual_spectrum_intensities, min_mz, max_mz, bin_size)

    norm_actual_vector = actual_spectrum_vector / np.linalg.norm(actual_spectrum_vector) if np.linalg.norm(actual_spectrum_vector) else actual_spectrum_vector
    norm_theoretical_vector = theoretical_spectrum_vector / np.linalg.norm(theoretical_spectrum_vector)

    score = cosine_similarity(norm_theoretical_vector, norm_actual_vector)
    return score


def match_spectrum_to_peptides(spectrum, peptides, min_mz, max_mz):
    best_match = None
    best_score = -1

    spectrum_pepmass = spectrum['params']['pepmass'][0]
    spectrum_charge = spectrum['params']['charge'][0]
    spectrum_mass = calculate_spectrum_mass(spectrum_pepmass, spectrum_charge)
    
    for peptide in peptides:
        if abs(spectrum_mass - peptides[peptide]) < 0.1: # optimization
            score = score_peptide_spectrum_match_intensity(peptide, spectrum, min_mz, max_mz)
            if score > best_score:
                best_score = score
                best_match = peptide
                
    return best_match, best_score

def calculate_fdr(matches, target_peptides):
    total_decoys = 0
    total_targets = 0
    largest_n = 0
    
    for i, match in enumerate(matches):
        if match[-1] in target_peptides:
            total_targets += 1
        else:
            total_decoys += 1
            
        if total_decoys / total_targets <= 0.01:
            largest_n = i + 1 
    
    return largest_n

def main():
    spectra = sys.argv[2]
    target = sys.argv[4]
    decoy = sys.argv[6]
    output = sys.argv[8]
    
    target_peptides = read_peptides(target)
    decoy_peptides = read_peptides(decoy)
    all_peptides = {**target_peptides, **decoy_peptides}
    # print(f"Number of target peptides: {len(target_peptides)}")
    # print(f"Number of decoy peptides: {len(decoy_peptides)}")

    matches = []
    with mgf.read(spectra) as result:
        spectra = list(result)
        global_min_mz = min([min(spectrum['m/z array']) for spectrum in spectra])
        global_max_mz = max([max(spectrum['m/z array']) for spectrum in spectra])
        
        print(f"Number of spectra: {len(spectra)}")
        for i, spectrum in enumerate(spectra):
            match, score = match_spectrum_to_peptides(spectrum, all_peptides, global_min_mz, global_max_mz)
            if match:
                matches.append((i, spectrum['params']['pepmass'][0], spectrum['params']['charge'][0], score, match))
    
    matches.sort(key=lambda x: x[3], reverse=True) # sort matches by score in descending order
    
    with open(output, 'w') as outfile:
        outfile.write('Id\tm/z\tz\tScore\tPeptide\n')
        for match in matches:
            outfile.write(f"{match[0]}\t{match[1]:.4f}\t{match[2]}\t{match[3]}\t{match[4]}\n")

    n = calculate_fdr(matches, target_peptides)
    print(f"Largest n: {n}")

if __name__ == "__main__":
    main()
