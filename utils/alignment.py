import subprocess, tempfile, os
from typing import List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def _call_muscle(input_fasta: str) -> str:
    with tempfile.TemporaryDirectory() as td:
        in_fa = os.path.join(td, "in.fasta")
        out_fa = os.path.join(td, "out.fasta")
        with open(in_fa, "w") as f:
            f.write(input_fasta)
        # MUSCLE v5 syntax
        cmd = ["muscle", "-align", in_fa, "-output", out_fa]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        with open(out_fa, "r") as f:
            return f.read()

def percent_identity(ref_aln: str, q_aln: str) -> float:
    # Identity over non-gap positions in the reference
    matches = 0
    denom = 0
    for r, q in zip(ref_aln, q_aln):
        if r != "-":
            denom += 1
            if q == r:
                matches += 1
    return 100.0 * matches / denom if denom else 0.0

def align_and_score(reference_seq: str, sequences: List[Tuple[str, str]]):
    """
    reference_seq: raw AA string for reference
    sequences: list of (name, raw AA string) for candidates
    Returns: list of dicts: name, length, pct_identity
    """
    # Build FASTA with reference first
    recs = [SeqRecord(Seq(reference_seq), id="REF", description="")]
    for name, seq in sequences:
        recs.append(SeqRecord(Seq(seq), id=name, description=""))
    with tempfile.TemporaryDirectory() as td:
        in_fa = os.path.join(td, "all.fasta")
        SeqIO.write(recs, in_fa, "fasta")
        with open(in_fa) as f:
            aln_out = _call_muscle(f.read())
    # Parse aligned FASTA
    aligned = list(SeqIO.parse(tempfile.NamedTemporaryFile(delete=False, mode="w+"), "fasta"))
