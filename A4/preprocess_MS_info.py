def parse_mgf(file_path):
    spectra = []
    with open(file_path, 'r') as file:
        spectrum = {}
        for line in file:
            line = line.strip()
            if line.startswith('BEGIN IONS'):
                spectrum = {}
            elif line.startswith('END IONS'):
                spectra.append(spectrum)
            elif line.startswith('PEPMASS='):
                spectrum['m/z'] = float(line.split('=')[1])
            elif line.startswith('CHARGE='):
                charge_info = line.split('=')[1].rstrip('+').rstrip('-')  # Remove + or - from the charge
                spectrum['z'] = int(charge_info)
            elif line.startswith('TITLE='):
                spectrum['title'] = line.split('=')[1]
    return spectra

# Example usage
file_path = './test1.mgf'
spectra_info = parse_mgf(file_path)
for spectrum in spectra_info:
    print(spectrum)
