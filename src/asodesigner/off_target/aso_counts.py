import subprocess
import random
from typing import Optional


def run_aso_counts(
        aso_fasta_file: str,
        k: int = 2,
        target_file: Optional[str] = None
) -> str:
    """
    Run aso_counts.sh and return only the string that appears after '__RETURN__:' in stdout.
    Raises RuntimeError if the script fails.

    Args:
        aso_fasta_file: Path to ASO query FASTA file
        k: Maximum Hamming distance (default: 2)
        target_file: Optional target FASTA file. If provided, runs on-target calculation.
                     If None, processes all chromosomes (off-target).

    Returns:
        Path to the output JSON file

    Raises:
        RuntimeError: If the bash script encounters an error
    """
    import os
    import time

    if not os.path.exists(aso_fasta_file):
        raise FileNotFoundError(f"ASO FASTA file does not exist: {aso_fasta_file}")

    # Check file size to ensure it's not empty
    file_size = os.path.getsize(aso_fasta_file)

    if file_size == 0:
        raise ValueError(f"ASO FASTA file is empty: {aso_fasta_file}")

    if target_file:
        if not os.path.exists(target_file):
            raise FileNotFoundError(f"Target file does not exist: {target_file}")

    # Small delay to ensure file system sync (sometimes needed on network drives)
    time.sleep(0.1)

    session_id = random.randint(1, 1_000_000)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    aso_counts_script = os.path.join(script_dir, "aso_counts.sh")

    # Build command
    cmd = [
        "bash",
        aso_counts_script,
        "-q", aso_fasta_file,
        "-k", str(k),
        "-s", str(session_id)
    ]

    # Add target file if provided
    if target_file:
        cmd.extend(["-t", target_file])

    result = subprocess.run(cmd, capture_output=True, text=True)

    # Parse output for errors, info, and return value
    error_msgs = []
    info_msgs = []
    return_value = None

    for line in result.stdout.splitlines():
        if line.startswith("__ERROR__:"):
            error_msgs.append(line.split(":", 1)[1].strip())
            print(f"[ERROR] {line}")
        elif line.startswith("__ERROR_LINE__:"):
            error_msgs.append(f"  Line: {line.split(':', 1)[1].strip()}")
            print(f"[ERROR] {line}")
        elif line.startswith("__ERROR_CMD__:"):
            error_msgs.append(f"  Command: {line.split(':', 1)[1].strip()}")
            print(f"[ERROR] {line}")
        elif line.startswith("__RETURN__:"):
            return_value = line.split(":", 1)[1].strip()
            print(f"[RETURN] {line}")
        else:
            # Print other output for debugging
            if line.strip():
                print(f"[BASH] {line}")

    # Print stderr if present
    if result.stderr:
        print(f"[STDERR] {result.stderr}")

    # Check if errors were found
    if error_msgs:
        error_text = "\n".join(error_msgs)
        raise RuntimeError(f"Bash script aso_counts.sh failed:\n{error_text}")

    # Check for non-zero return code
    if result.returncode != 0:
        raise RuntimeError(
            f"Bash script aso_counts.sh exited with code {result.returncode}\n"
            f"STDOUT: {result.stdout}\n"
            f"STDERR: {result.stderr}"
        )

    # Check if return value was found
    if not return_value:
        raise RuntimeError(
            f"No __RETURN__ marker found in bash script output.\n"
            f"This usually means the script completed but didn't produce the expected output.\n"
            f"Last output: {result.stdout[-500:] if result.stdout else '(empty)'}"
        )

    print(f"[DEBUG run_aso_counts] Success! Output file: {return_value}")
    return return_value


# Example usage:
if __name__ == "__main__":
    try:
        # Off-target (all chromosomes)
        print("=" * 60)
        print("TEST 1: Off-target (all chromosomes)")
        print("=" * 60)
        output_path = run_aso_counts("aso_query.fa")
        print(f"✓ Off-target results: {output_path}")
        print()

        # On-target (specific gene)
        print("=" * 60)
        print("TEST 2: On-target (specific gene)")
        print("=" * 60)
        output_path = run_aso_counts("aso_query.fa", k=2, target_file="target_gene.fa")
        print(f"✓ On-target results: {output_path}")

    except RuntimeError as e:
        print(f"✗ Error occurred:")
        print(str(e))
    except Exception as e:
        print(f"✗ Unexpected error:")
        print(f"{type(e).__name__}: {e}")
