from pyteomics import mgf, mass
import numpy as np
import sys

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
        b_ion = mass.fast_mass(peptide[:i], ion_type='b')
        y_ion = mass.fast_mass(peptide[i:], ion_type='y')
        ions.extend([b_ion, y_ion])
    return np.array(sorted(ions))

def cosine_similarity(v1, v2):
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    return dot_product / (norm_v1 * norm_v2) if norm_v1 and norm_v2 else 0

# check if mass difference is within tolerance
# if yes, evaluate score
def score_peptide_spectrum_match_intensity(peptide, spectrum):
    theoretical_spectrum = generate_theoretical_spectrum(peptide)
    actual_spectrum_mz = np.array(spectrum['m/z array'])
    actual_spectrum_intensities = np.array(spectrum['intensity array'])
    
    # Initialize a vector for the actual spectrum based on theoretical spectrum size
    actual_spectrum_vector = np.zeros_like(theoretical_spectrum)
    tolerance = 0.5  # Mass tolerance for matching peaks
    
    # Assign intensities to the actual spectrum vector based on matching m/z values
    for mz, intensity in zip(actual_spectrum_mz, actual_spectrum_intensities):
        diff = np.abs(theoretical_spectrum - mz)
        closest_index = np.argmin(diff)
        if diff[closest_index] <= tolerance:
            actual_spectrum_vector[closest_index] += intensity  # Sum intensities for simplicity
    
    # Normalize vectors to have a unit norm before calculating cosine similarity
    norm_actual_vector = actual_spectrum_vector / np.linalg.norm(actual_spectrum_vector) if np.linalg.norm(actual_spectrum_vector) else actual_spectrum_vector
    theoretical_spectrum_vector = np.ones_like(theoretical_spectrum)  # Assume uniform intensity for theoretical peaks
    norm_theoretical_vector = theoretical_spectrum_vector / np.linalg.norm(theoretical_spectrum_vector)
    
    # Calculate cosine similarity with intensity information
    score = cosine_similarity(norm_theoretical_vector, norm_actual_vector)
    return score


def match_spectrum_to_peptides(spectrum, peptides):
    best_match = None
    best_score = -1  # Initialize with a score that any match will exceed

    spectrum_pepmass = spectrum['params']['pepmass'][0]
    spectrum_charge = spectrum['params']['charge'][0]
    spectrum_mass = calculate_spectrum_mass(spectrum_pepmass, spectrum_charge)
    
    for peptide in peptides:
        if abs(spectrum_mass - peptides[peptide]) < 0.1:
            score = score_peptide_spectrum_match_intensity(peptide, spectrum)
            if score > best_score:
                best_score = score
                best_match = peptide
                
    return best_match, best_score

def calculate_fdr(matches):
    total_decoys = 0
    total_targets = 0
    largest_n = 0
    
    for i, match in enumerate(matches):
        if match[-1] == 'decoy':
            total_decoys += 1
        else:
            total_targets += 1
        
        # Calculate FDR up to the current position
        if total_targets > 0:
            current_fdr = total_decoys / (total_targets + total_decoys)
        else:
            current_fdr = 0
        
        # Check if current FDR is within the acceptable limit (<= 1%)
        if current_fdr <= 0.01:
            largest_n = i + 1  # Adjust for 0-based indexing
    
    return largest_n

def main():
    spectra = sys.argv[2]
    target = sys.argv[4]
    decoy = sys.argv[6]
    output = sys.argv[8]
    
    # spectra = parse_mgf("./test1.mgf")
    target_peptides = read_peptides(target)
    decoy_peptides = read_peptides(decoy)
    all_peptides = {**target_peptides, **decoy_peptides}
    # print(f"Number of target peptides: {len(target_peptides)}")
    # print(f"Number of decoy peptides: {len(decoy_peptides)}")

    matches = []
    with mgf.read(spectra) as spectra:
        print(f"Number of spectra: {len(spectra)}")
        for i, spectrum in enumerate(spectra):
            match, score = match_spectrum_to_peptides(spectrum, all_peptides)
            if match:
                matches.append((i, spectrum['params']['pepmass'][0], spectrum['params']['charge'][0], score, match))
    
    matches.sort(key=lambda x: x[3], reverse=True) # sort matches by score in descending order
    
    with open(output, 'w') as outfile:
        outfile.write('Id\tm/z\tz\tScore\tPeptide\n')
        for match in matches:
            outfile.write(f"{match[0]}\t{match[1]:.4f}\t{match[2]}\t{match[3]}\t{match[4]}\n")

    n = calculate_fdr(matches)
    print(f"Largest n: {n}")

if __name__ == "__main__":
    main()
