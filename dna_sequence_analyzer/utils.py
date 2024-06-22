
# Define any utilities that may be needed

def parse_sequences(filename):
    """
    Parses a FASTA file and extracts DNA sequences into a dictionary.

    Args:
        - filename (str): Path to the FASTA file containing DNA sequences.
    
    Returns:
        - dict: A dictionary where keys are sequence IDs (from FASTA headers)
          and values are corresponding DNA sequences.
    """
    sequences = {}
    current_seq_id = ''
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_seq_id = line[1:]
                sequences[current_seq_id] = ''
            else:
                sequences[current_seq_id] += line
    return sequences