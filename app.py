import gradio as gr
import pandas as pd
from protein_mpnn_wrapper import run_protein_mpnn  # This is our mock script
from utils.alignment import run_muscle, compute_percent_identity, alignment_to_clustal_text

def process_and_align(pdb_file, chains_to_design, num_sequences):
    """
    Main function for the Gradio interface.
    1. Runs ProteinMPNN (mocked).
    2. Aligns the results with the original sequence.
    3. Computes homology.
    4. Returns results for display.
    """
    if pdb_file is None or not chains_to_design:
        return "Please upload a PDB file and specify chains.", None, None, None

    try:
        # 1. Run ProteinMPNN to get designed sequences and the reference sequence
        # Note: The wrapper needs to be adapted to parse the PDB for the native sequence.
        # For this example, we'll use a placeholder native sequence.
        native_sequence = "PIAQIHILEGRSDEQKETLIREVSEAISRSLDAPLTSVRVIITEMAKGHFGIGGELASK" # Placeholder
        designed_seqs_dict = run_protein_mpnn(pdb_file.name, chains_to_design, num_sequences)

        if not designed_seqs_dict:
            return "ProteinMPNN returned no sequences.", None, None, None

        # 2. Prepare sequences for alignment
        fasta_data = {"native_reference": native_sequence}
        candidate_seqs = []
        for i, seq in enumerate(designed_seqs_dict.values()):
            seq_id = f"designed_seq_{i+1}"
            fasta_data[seq_id] = seq
            candidate_seqs.append((seq_id, seq))

        # 3. Run MUSCLE alignment
        alignment = run_muscle(fasta_data)
        alignment_text = alignment_to_clustal_text(alignment)

        # 4. Compute percent identity
        homology_results = compute_percent_identity(alignment, reference_id="native_reference")

        # 5. Format results for Gradio output
        homology_df = pd.DataFrame(
            list(homology_results.items()),
            columns=["Sequence ID", "Percent Identity vs Native"]
        ).sort_values(by="Percent Identity vs Native", ascending=False)
        
        designed_seqs_text = "\n".join(f">{key}\n{value}" for key, value in designed_seqs_dict.items())

        return designed_seqs_text, alignment_text, homology_df, "Processing complete."

    except Exception as e:
        return None, None, None, f"An error occurred: {str(e)}"

# --- Gradio Interface ---
with gr.Blocks() as demo:
    gr.Markdown("# ProteinMPNN with Sequence Alignment")
    gr.Markdown("Design protein sequences with ProteinMPNN and see how they align against the native sequence.")

    with gr.Row():
        with gr.Column(scale=1):
            pdb_file = gr.File(label="Upload PDB File")
            chains_to_design = gr.Textbox(label="Chains to Design", placeholder="e.g., A,B")
            num_sequences = gr.Slider(label="Number of Sequences to Generate", minimum=1, maximum=20, step=1, value=4)
            run_button = gr.Button("Run Design and Alignment", variant="primary")

        with gr.Column(scale=2):
            status_output = gr.Textbox(label="Status", interactive=False)
            designed_seqs_output = gr.Textbox(label="Designed Sequences (FASTA)", lines=8, interactive=False)

    with gr.Row():
        with gr.Accordion("Alignment and Homology Results", open=True):
            alignment_output = gr.Textbox(label="Sequence Alignment (MUSCLE)", lines=15, interactive=False)
            homology_output = gr.DataFrame(label="Percent Identity vs. Native", interactive=False)
    
    run_button.click(
        fn=process_and_align,
        inputs=[pdb_file, chains_to_design, num_sequences],
        outputs=[designed_seqs_output, alignment_output, homology_output, status_output]
    )

if __name__ == "__main__":
    demo.launch()