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
from backend_utils import file_to_id
from aso_generator import foo
from send_email import send_processing_started_email, send_processing_completed_email, prepare_data_before_sending

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
        gene_data = request_data.geneSequence
        organism_name = file_to_id[request_data.organismFile]
        
        logger.info(f"Running ASO generation for gene: {request_data.geneName}, organism: {organism_name}")
        
        # foo returns: (dataframe, full_mrna_sequence, gene_name, organism_name)
        # or: (dataframe, full_mrna_sequence, gene_name)
        # or: (dataframe, full_mrna_sequence)
        result = foo(
            organism_name,
            request_data.geneName,
            gene_data,
            request_data.top_k,
            request_data.includeFeatureBreakdown,
        )
        
        logger.info(f"ASO generation completed. Result type: {type(result)}, length: {len(result) if hasattr(result, '__len__') else 'N/A'}")
        
        # Ensure result is a tuple with at least (df, full_mrna_seq)
        # Augment with gene_name and organism_name if not provided by foo
        if isinstance(result, tuple):
            if len(result) == 2:
                # Add gene_name and organism_name
                result = (*result, request_data.geneName, organism_name)
            elif len(result) == 3:
                # Add organism_name
                result = (*result, organism_name)
        else:
            # If foo returns just a dataframe, wrap it
            logger.warning("foo returned unexpected format, attempting to handle gracefully")
            result = (result, "", request_data.geneName, organism_name)
        
        # Prepare comprehensive files for biologists
        files = prepare_data_before_sending(result, website_url="https://2025.igem.wiki/tau-israel/model")
        
        logger.info(f"Prepared {len(files)} files for email delivery")

        # Send "processing completed" email with all attachments
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
            files=files
        ))
        
        logger.info(f"Successfully sent completion email to {request_data.userEmail}")

    except Exception as e:
        logger.exception("Error in background ASO task")
        # Optionally, send an error email to the user
        try:
            asyncio.run(send_error_email(
                to=request_data.userEmail,
                name=request_data.userName,
                error_message=str(e)
            ))
        except:
            logger.error("Failed to send error notification email")


async def send_error_email(to: str, name: str, error_message: str):
    """Send email notification when processing fails."""
    from send_email import send_mail, escape_html
    
    subject = "ASO Analysis Error"
    html = f"""
    <div style="font-family: Arial, sans-serif; line-height: 1.6;">
      <h2>Hello {escape_html(name)},</h2>
      <p>We encountered an error while processing your ASO analysis request.</p>
      <p><strong>Error details:</strong></p>
      <pre style="background: #f5f5f5; padding: 10px; border-radius: 5px; overflow-x: auto;">{escape_html(error_message)}</pre>
      <p>Please check your input parameters and try again. If the problem persists, contact support.</p>
      <p style="color:#555; font-size: 12px;">This message was sent from Oncoligo ASO Generator.</p>
    </div>
    """
    await send_mail(to, subject, html)


@app.post("/run_aso")
async def run_aso(request_data: RunASORequest, background_tasks: BackgroundTasks):
    """
    Main endpoint to run ASO generation.
    Sends immediate confirmation email and processes ASO generation in background.
    """
    logger.info(f"Received ASO request from {request_data.userEmail} for gene {request_data.geneName}")
    
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
    logger.info(f"ASO processing started in background for {request_data.userEmail}")
    return {"status": "Processing started", "message": "You will receive an email when the analysis is complete"}


@app.get("/")
async def root():
    """Health check endpoint."""
    return {"status": "ok", "message": "ASO Generator API"}


@app.get("/health")
async def health_check():
    """Detailed health check endpoint."""
    return {
        "status": "healthy",
        "service": "Oncoligo ASO Generator",
        "version": "1.0.0"
    }