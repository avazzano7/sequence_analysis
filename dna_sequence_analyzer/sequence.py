
import sys
import os

# Add the parent directory of 'dna_sequence_analyzer' to sys.path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

# Now you can import parse_sequences from dna_sequence_analyzer.utils
from dna_sequence_analyzer.utils import parse_sequences

def calculate_gc_content(sequence):
    """
    Calculate the GC content percentage of a DNA sequence.
    """
    total_bases = len(sequence)
    sequence = sequence.lower()
    gc_count = sequence.count('g') + sequence.count('c')
    if total_bases > 0:
        gc_content = (gc_count / total_bases) * 100
    else:
        gc_content = 0.0
    return gc_content


def translate_to_protein(sequence):
    """
    Translate a DNA sequence into a protein sequence using a codon table.
    """
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    protein_sequence = ''
    for i in range(0, len(sequence) - 2, 3):
        codon = sequence[i:i+3].upper()  # Ensure codon is in uppercase
        if codon in codon_table:
            protein_sequence += codon_table[codon]
        else:
            protein_sequence += 'X'  # Unknown codon represented as 'X'
    return protein_sequence


def analyze_sequence(filename):
    """
    Analyze DNA sequences from a FASTA file and print sequence details.
    """
    try:
        sequences = parse_sequences(filename)
        for seq_id, sequence in sequences.items():
            sequence = sequence.replace('\n', '')  # Remove newline characters if any
            length = len(sequence)
            gc_content = calculate_gc_content(sequence)

            print(f"Sequence ID: {seq_id}")
            print(f"Sequence length: {length} nucleotides")
            print(f"GC content: {gc_content:.2f}%")

            # Example translation to protein sequence
            protein_sequence = translate_to_protein(sequence)
            print(f"\nProtein sequence:\n{protein_sequence}")

    except Exception as e:
        print(f"Error: {e}")


def main():
    import sys
    if len(sys.argv) != 2:
        print("Usage: dna-analyzer <sequence_file>")
        sys.exit(1)
    filename = sys.argv[1]
    analyze_sequence(filename)

if __name__ == '__main__':
    main()
