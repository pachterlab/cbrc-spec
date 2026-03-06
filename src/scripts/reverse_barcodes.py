bases = {"A": "T", "T": "A", "C": "G", "G": "C"}

def reverse_barcode(barcode):
    """Return the reverse complement of a DNA barcode string.

    Args:
        barcode (str): A string representing a DNA barcode (e.g., "ATCG").

    Returns:
        str: The reverse complement of the input barcode.
    """
    reversed_barcode = barcode[::-1]
    reverse_complement = ''.join(bases[base] for base in reversed_barcode)
    return reverse_complement

path = "seqspec/onlists/Parse"
files = [f"{path}/v2/BC1.txt",
         f"{path}/v2/BC2_3.txt",
         f"{path}/v3/BC1.txt",
         f"{path}/v3/BC2_3.txt"]

for file in files:
    with open(file, 'r') as infile, open(file.replace('.txt', '_rev.txt'), 'w') as outfile:
        for line in infile:
            barcode = line.strip()
            revcomp_barcode = reverse_barcode(barcode)
            outfile.write(revcomp_barcode + '\n')