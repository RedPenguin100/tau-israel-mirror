import subprocess
import random
from typing import Optional
import os

def run_aso_counts(
    aso_fasta_file: str,
    k: int = 2,
    target_file: Optional[str] = None
) -> str:
    """
    Run aso_counts.sh and return only the string that appears after '__RETURN__:' in stdout.
    
    Args:
        aso_fasta_file: Path to ASO query FASTA file
        k: Maximum Hamming distance (default: 2)
        target_file: Optional target FASTA file. If provided, runs on-target calculation.
                     If None, processes all chromosomes (off-target).
    
    Returns:
        Path to the output JSON file
    """
    session_id = random.randint(1, 1_000_000)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    aso_counts_script = os.path.join(script_dir, "aso_counts.sh")
    
    # Build command
    cmd = [
        "bash", 
        "aso_counts.sh", 
        "-q", aso_fasta_file, 
        "-k", str(k), 
        "-s", str(session_id)
    ]
    
    # Add target file if provided
    if target_file:
        cmd.extend(["-t", target_file])
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=script_dir)
    
    # Parse output for return value
    for line in result.stdout.splitlines():
        if line.startswith("__RETURN__:"):
            return line.split(":", 1)[1].strip()
    
    return ""  # if no marker found


# Example usage:
if __name__ == "__main__":
    # Off-target (all chromosomes)
    output_path = run_aso_counts("aso_query.fa")
    print(f"Off-target results: {output_path}")
    
    # On-target (specific gene)
    output_path = run_aso_counts("aso_query.fa", k=2, target_file="target_gene.fa")
    print(f"On-target results: {output_path}")
