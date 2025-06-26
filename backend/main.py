from fastapi import FastAPI, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from datetime import datetime
from pathlib import Path
from Bio import SeqIO
import shutil
from fastapi import FastAPI, Request
from pydantic import BaseModel
from typing import Optional, Dict, Any
from fastapi.middleware.cors import CORSMiddleware
import os

app = FastAPI()

# CORS setup
ORIGINS = [
    "http://localhost:3000",
    "https://2025.igem.wiki",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

UPLOAD_FOLDER = Path("uploads")
UPLOAD_FOLDER.mkdir(exist_ok=True)



import os
import time

def remove_expired_fasta_files(folder_path=UPLOAD_FOLDER, max_age_ms=18640):
    """
    Remove .fasta files with names in format: {file.name}_{timestamp}.fasta
    if the timestamp is older than max_age_ms (milliseconds).
    
    Args:
        folder_path (str): Path to the folder containing the files.
        max_age_ms (int): Maximum allowed age in milliseconds.
    """
    now = int(time.time() * 1000)

    for filename in os.listdir(folder_path):
        if filename.endswith('.fasta'):
            parts = filename.rsplit('_', 1)
            if len(parts) == 2 and parts[1].endswith('.fasta'):
                timestamp_str = parts[1][:-6]  # remove '.fasta'
                if timestamp_str.isdigit():
                    file_time = int(timestamp_str)
                    age = now - file_time
                    if age > max_age_ms:
                        file_path = os.path.join(folder_path, filename)
                        try:
                            os.remove(file_path)
                            print(f"Deleted: {file_path}")
                        except Exception as e:
                            print(f"Failed to delete {file_path}: {e}")


@app.get("/")
def index():
    return {
        "message": "Welcome to the Gene Upload API! Use /upload_gene_fa or /upload_gene_gb to upload files in chunks."
    }

async def save_chunk(filename: Path, chunk_file: UploadFile):
    """Append chunk content to file."""
    mode = "ab"  # append in binary mode
    with open(filename, mode) as f:
        shutil.copyfileobj(chunk_file.file, f)
    await chunk_file.close()

def validate_fasta(file_path: Path):
    try:
        list(SeqIO.parse(str(file_path), "fasta"))
    except Exception as e:
        raise ValueError(f"Invalid FASTA format: {str(e)}")

def validate_genbank(file_path: Path):
    try:
        list(SeqIO.parse(str(file_path), "genbank"))
    except Exception as e:
        raise ValueError(f"Invalid GenBank format: {str(e)}")

@app.post("/upload_gene_fa")
async def upload_gene_fa(
    filename: str = Form(...),
    is_last_chunk: str = Form(...),
    chunk: UploadFile = File(...)
):
    try:
        save_path = UPLOAD_FOLDER / filename

        # Save chunk (append mode)
        await save_chunk(save_path, chunk)

        # If this is the last chunk, validate file
        if is_last_chunk.lower() == "true":
            try:
                remove_expired_fasta_files()
                validate_fasta(save_path)
            except ValueError as ve:
                # Delete invalid file to avoid confusion
                save_path.unlink(missing_ok=True)
                return JSONResponse(
                    status_code=400,
                    content={"error": str(ve)}
                )
        return {"message": "Chunk saved"}

    except Exception as e:
        return JSONResponse(status_code=500, content={"error": str(e)})


@app.post("/upload_gene_gb")
async def upload_gene_gb(
    filename: str = Form(...),
    is_last_chunk: str = Form(...),
    chunk: UploadFile = File(...)
):
    try:
        save_path = UPLOAD_FOLDER / filename

        # Save chunk (append mode)
        await save_chunk(save_path, chunk)

        # If this is the last chunk, validate file
        if is_last_chunk.lower() == "true":
            try:
                remove_expired_fasta_files()
                validate_genbank(save_path)
            except ValueError as ve:
                save_path.unlink(missing_ok=True)
                return JSONResponse(
                    status_code=400,
                    content={"error": str(ve)}
                )
        return {"message": "Chunk saved"}

    except Exception as e:
        return JSONResponse(status_code=500, content={"error": str(e)})
# Pydantic model for input validation
class GenASORequest(BaseModel):
    geneInputType: str
    organismFile: str
    geneFile: Optional[str] = None
    geneSequence: Optional[str] = None
    numericParams: Dict[str, Any]
    viewASO: bool

# Placeholder outer function
def foo(organismFile, geneData, numericParams, viewASO):
    """
    geneData is either geneFile (str) or geneSequence (str) depending on geneInputType.
    This function returns a dict with 'asoSequence': list of {name, sequence}.
    """
    # Dummy implementation for testing
    return {
        "asoSequence": [
            {"name": "ASO_1", "sequence": "ACGTACGTACGT"},
            {"name": "ASO_2", "sequence": "TTGGAACC"}
        ]
    }

@app.post("/gen_aso")
async def generate_aso(request_data: GenASORequest):
    try:
        if request_data.geneInputType == "Upload Gene Fasta":
            if not request_data.geneFile:
                return {"error": "geneFile is required for 'Upload Gene Fasta'"}
            result = foo(
                request_data.organismFile,
                request_data.geneFile,
                request_data.numericParams,
                request_data.viewASO,
            )
        elif request_data.geneInputType == "Enter Gene Sequence":
            if not request_data.geneSequence:
                return {"error": "geneSequence is required for 'Enter Gene Sequence'"}
            result = foo(
                request_data.organismFile,
                request_data.geneSequence,
                request_data.numericParams,
                request_data.viewASO,
            )
        else:
            return {"error": f"Unsupported geneInputType: {request_data.geneInputType}"}

        if "asoSequence" not in result:
            return {"error": "asoSequence not returned from processing"}

        return result

    except Exception as e:
        return {"error": f"Internal server error: {str(e)}"}