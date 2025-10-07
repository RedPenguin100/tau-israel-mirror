"""
Main FastAPI Application
Handles HTTP endpoints and orchestrates ASO generation workflow.
"""
import base64
import logging
import asyncio
from typing import Optional
from fastapi import FastAPI, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO

from aso_gen import foo
from send_email import send_processing_started_email, send_processing_completed_email

# Initialize FastAPI app
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

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# Request models
class RunASORequest(BaseModel):
    geneName: str
    organismFile: str
    geneFile: Optional[str] = None
    geneSequence: Optional[str] = None
    top_k: int
    includeFeatureBreakdown: bool
    userEmail: str
    userName: str


def run_aso_task(request_data: RunASORequest):
    """Background task to run ASO generation and send completion email."""
    try:
        # Generate ASO sequences (long-running)
        gene_data = request_data.geneSequence or request_data.geneFile or request_data.geneName
        result = foo(
            request_data.organismFile,
            request_data.geneName,
            gene_data,
            request_data.top_k,
            request_data.includeFeatureBreakdown,
        )

        # Convert ASO to FASTA
        fasta_records = [
            SeqRecord(Seq(aso["sequence"]), id=aso["name"], description="") 
            for aso in result["asoSequence"]
        ]
        fasta_str_io = StringIO()
        SeqIO.write(fasta_records, fasta_str_io, "fasta")
        fasta_bytes = fasta_str_io.getvalue().encode()
        fasta_base64 = base64.b64encode(fasta_bytes).decode()

        # Send "processing completed" email
        asyncio.run(send_processing_completed_email(
            to=request_data.userEmail,
            name=request_data.userName,
            asoData={
                "geneName": request_data.geneName,
                "organismFile": request_data.organismFile,
                "geneFile": request_data.geneFile,
                "geneSequence": request_data.geneSequence,
                "top_k": request_data.top_k,
                "includeFeatureBreakdown": request_data.includeFeatureBreakdown,
            },
            file=fasta_base64
        ))

    except Exception as e:
        logging.exception("Error in background ASO task")


@app.post("/run_aso")
async def run_aso(request_data: RunASORequest, background_tasks: BackgroundTasks):
    """
    Main endpoint to run ASO generation.
    Sends immediate confirmation email and processes ASO generation in background.
    """
    # Send "processing started" email immediately
    await send_processing_started_email(
        to=request_data.userEmail,
        name=request_data.userName,
        asoData={
            "geneName": request_data.geneName,
            "organismFile": request_data.organismFile,
            "geneFile": request_data.geneFile,
            "geneSequence": request_data.geneSequence,
            "top_k": request_data.top_k,
            "includeFeatureBreakdown": request_data.includeFeatureBreakdown,
        },
    )

    # Add the long-running ASO generation to background tasks
    background_tasks.add_task(run_aso_task, request_data)

    # Return immediately
    return {"status": "Processing started"}


@app.get("/")
async def root():
    """Health check endpoint."""
    return {"status": "ok", "message": "ASO Generator API"}
