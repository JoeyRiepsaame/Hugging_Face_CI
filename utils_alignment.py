import subprocess
from io import StringIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from typing import Dict

def run_muscle(fasta_dict: Dict[str, str]) -> MultipleSeqAlignment:
    """Runs MUSCLE on a dictionary of sequences."""
    fasta_str = ""
    for seq_id, sequence in fasta_dict.items():
        fasta_str += f">{seq_id}\n{sequence}\n"
    
    # Run MUSCLE process
    process = subprocess.run(
        ['muscle', '-in', '-', '-out', '-'],
        input=fasta_str.encode('utf-8'),
        capture_output=True,
        check=True
    )
    
    # Parse the aligned output
    alignment = AlignIO.read(StringIO(process.stdout.decode('utf-8')), "fasta")
    return alignment

def compute_percent_identity(alignment: MultipleSeqAlignment, reference_id: str) -> Dict[str, float]:
    """Computes percent identity of each sequence against a reference sequence."""
    identities = {}
    try:
        ref_record = next(rec for rec in alignment if rec.id == reference_id)
    except StopIteration:
        raise ValueError(f"Reference ID '{reference_id}' not found in alignment.")

    ref_seq = str(ref_record.seq)
    
    for record in alignment:
        if record.id == reference_id:
            continue
        
        matches = 0
        total_len = 0
        comp_seq = str(record.seq)

        for i in range(len(ref_seq)):
            # Only compare positions where the reference is not a gap
            if ref_seq[i] != '-':
                total_len += 1
                if ref_seq[i] == comp_seq[i]:
                    matches += 1
        
        percent_identity = (matches / total_len) * 100 if total_len > 0 else 0
        identities[record.id] = round(percent_identity, 2)
        
    return identities

def alignment_to_clustal_text(alignment: MultipleSeqAlignment) -> str:
    """Formats the alignment into a string similar to CLUSTAL format."""
    return alignment.format("clustal")