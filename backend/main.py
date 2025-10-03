import base64
import logging
import os
import time
import uuid
from datetime import datetime
from email.message import EmailMessage
from pathlib import Path
from typing import Optional

from Bio import SeqIO
from dotenv import load_dotenv
from fastapi import FastAPI, UploadFile, File, Form, Request, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel
import aiosmtplib


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

BASE_DIR = Path(__file__).resolve().parent.parent

UPLOAD_FOLDER = Path("uploads")
UPLOAD_FOLDER.mkdir(exist_ok=True)

GENE_NAME_CANDIDATES = [
    BASE_DIR / "client" / "public" / "gene_names.txt",
    BASE_DIR / "tests" / "gene_names.txt",
]


def _resolve_gene_names_path() -> Optional[Path]:
    for candidate in GENE_NAME_CANDIDATES:
        if candidate.exists():
            return candidate
    return None


class GenASORequest(BaseModel):
    geneName: str
    organismFile: str
    geneFile: Optional[str] = None
    geneSequence: Optional[str] = None
    top_k: int
    includeFeatureBreakdown: bool


# Placeholder outer function
def foo(organismFile, geneName, geneData, top_k, includeFeatureBreakdown):
    """
    geneName: selected gene identifier.
    geneData: nucleotide sequence or identifier supplied by the client (optional).
    This function returns a dict with 'asoSequence': list of {name, sequence}.
    """
    # Dummy implementation for testing
    time.sleep(30)
    return {
        "asoSequence": [
            {"name": "ASO_1", "sequence": "ACGTACGTACGT"},
            {"name": "ASO_2", "sequence": "TTGGAACC"}
        ]
    }


@app.post("/gen_aso")
async def generate_aso(request_data: GenASORequest):
    try:
        gene_data = request_data.geneSequence or request_data.geneFile or request_data.geneName

        result = foo(
            request_data.organismFile,
            request_data.geneName,
            gene_data,
            request_data.top_k,
            request_data.includeFeatureBreakdown,
        )

        if "asoSequence" not in result:
            return {"error": "asoSequence not returned from processing"}

        return result

    except Exception as e:
        return {"error": f"Internal server error: {str(e)}"}


@app.get("/gene_names")
async def get_gene_names():
    try:
        gene_path = _resolve_gene_names_path()
        if not gene_path:
            raise HTTPException(status_code=500, detail="Gene names file not found")

        with gene_path.open("r", encoding="utf-8") as handle:
            names = [line.strip() for line in handle if line.strip()]

        return {"genes": names}
    except HTTPException:
        raise
    except Exception:
        logging.exception("Failed to read gene names")
        raise HTTPException(status_code=500, detail="Unable to load gene names")


# Load environment variables
load_dotenv(".env.local")

GMAIL_USER = os.getenv("GMAIL_USER")
GMAIL_PASS = os.getenv("GMAIL_APP_PASSWORD", "")
FROM_HEADER = os.getenv("EMAIL_FROM", f"Oncoligo ASO <{GMAIL_USER}>")

if not GMAIL_USER:
    raise RuntimeError("GMAIL_USER not set in environment variables")


class EmailRequest(BaseModel):
    email: str
    name: str
    type: str
    asoData: dict | None = {}
    file: str | None = None   # base64-encoded string


@app.post("/send-email")
async def send_email(req: EmailRequest):
    if not req.email or not req.name or not req.type:
        raise HTTPException(status_code=400, detail="Missing required fields")

    if not GMAIL_PASS:
        raise HTTPException(status_code=500, detail="GMAIL_APP_PASSWORD not set on server")

    try:
        if req.type == "processing_started":
            await send_processing_started_email(req.email, req.name, req.asoData or {})
        elif req.type == "processing_completed":
            await send_processing_completed_email(req.email, req.name, req.asoData or {}, req.file)
        else:
            raise HTTPException(status_code=400, detail="Unknown email type")

        return JSONResponse({"success": True})

    except Exception as e:
        logging.exception("Error sending email")
        raise HTTPException(status_code=500, detail="Failed to send email")


async def send_processing_started_email(to: str, name: str, asoData: dict):
    subject = "Your ASO analysis is being processed"

    organism_file = asoData.get("organismFile", "N/A")
    gene_name = asoData.get("geneName", "N/A")
    gene_seq = asoData.get("geneSequence") or ""
    gene_preview = gene_seq[:60] + ("…" if len(gene_seq) > 60 else "")
    top_k = asoData.get("top_k", "N/A")
    include_breakdown = "Yes" if asoData.get("includeFeatureBreakdown") else "No"

    html = f"""
    <div style="font-family: Arial, sans-serif; line-height: 1.6;">
      <h2>Hello {escape_html(name)},</h2>
      <p>We've received your ASO analysis request and processing has begun.</p>
      <h3>Summary of your request</h3>
      <ul>
        <li><strong>Organism file</strong>: {escape_html(str(organism_file))}</li>
        <li><strong>Gene name</strong>: {escape_html(str(gene_name))}</li>
        <li><strong>Gene sequence preview</strong>: <code>{escape_html(gene_preview or "N/A")}</code></li>
        <li><strong>Top-k results</strong>: {escape_html(str(top_k))}</li>
        <li><strong>Detailed analysis</strong>: {escape_html(include_breakdown)}</li>
      </ul>
      <p>You'll receive another email with your results once the analysis is complete.</p>
      <p style="color:#555; font-size: 12px;">This message was sent from Oncoligo ASO Generator.</p>
    </div>
    """

    await send_mail(to, subject, html)


async def send_processing_completed_email(to: str, name: str, asoData: dict, file: str | None):
    subject = "Your ASO analysis is complete!"

    organism_file = asoData.get("organismFile", "N/A")
    gene_name = asoData.get("geneName", "N/A")
    gene_seq = asoData.get("geneSequence") or ""
    gene_preview = gene_seq[:60] + ("…" if len(gene_seq) > 60 else "")
    top_k = asoData.get("top_k", "N/A")
    include_breakdown = "Yes" if asoData.get("includeFeatureBreakdown") else "No"

    html = f"""
    <div style="font-family: Arial, sans-serif; line-height: 1.6;">
      <h2>Hello {escape_html(name)},</h2>
      <p>Your ASO analysis has been completed successfully!</p>
      <h3>Summary</h3>
      <ul>
        <li><strong>Organism file</strong>: {escape_html(str(organism_file))}</li>
        <li><strong>Gene name</strong>: {escape_html(str(gene_name))}</li>
        <li><strong>Gene sequence preview</strong>: <code>{escape_html(gene_preview or "N/A")}</code></li>
        <li><strong>Top-k results</strong>: {escape_html(str(top_k))}</li>
        <li><strong>Detailed analysis</strong>: {escape_html(include_breakdown)}</li>
      </ul>
      <p style="color:#555; font-size: 12px;">This message was sent from Oncoligo ASO Generator.</p>
    </div>
    """

    attachments = []
    if file:
        try:
            decoded = base64.b64decode(file)
            attachments.append(("aso_results.fasta", decoded))
        except Exception:
            logging.warning("Failed to decode file, skipping attachment")

    await send_mail(to, subject, html, attachments)

import random

INVISIBLE_CHARS = ["\u200b", "\u200c", "\u200d", "\u2060", "\ufeff"]

def get_random_invisible(length: int = 3) -> str:
    """Return a random invisible string of given length."""
    return "".join(random.choice(INVISIBLE_CHARS) for _ in range(length))


async def send_mail(to: str, subject: str, html: str, attachments: list[tuple[str, bytes]] = []):
    msg = EmailMessage()
    msg["From"] = FROM_HEADER
    msg["To"] = to
    msg["Subject"] = subject + get_random_invisible()  # zero-width space to make Gmail treat it as new
    msg["Message-ID"] = f"<{uuid.uuid4()}@oncoligo.com>"
    msg.set_content("This email requires an HTML-capable client.")
    msg.add_alternative(html, subtype="html")

    for filename, content in attachments:
        msg.add_attachment(content, maintype="application", subtype="octet-stream", filename=filename)

    await aiosmtplib.send(
        msg,
        hostname="smtp.gmail.com",
        port=465,
        use_tls=True,
        username=GMAIL_USER,
        password=GMAIL_PASS,
    )


def escape_html(s: str) -> str:
    return (
        str(s)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&#039;")
    )
from fastapi import BackgroundTasks
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO


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
        fasta_records = [SeqRecord(Seq(aso["sequence"]), id=aso["name"], description="") for aso in result["asoSequence"]]
        fasta_str_io = StringIO()
        SeqIO.write(fasta_records, fasta_str_io, "fasta")
        fasta_bytes = fasta_str_io.getvalue().encode()
        fasta_base64 = base64.b64encode(fasta_bytes).decode()

        # Send "processing completed" email
        import asyncio
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
