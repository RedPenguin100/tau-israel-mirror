import os
import base64
import io
import zipfile
import logging
import random
from dotenv import load_dotenv
from brevo_python import Configuration, ApiClient, TransactionalEmailsApi
from brevo_python import SendSmtpEmail, SendSmtpEmailSender, SendSmtpEmailTo, SendSmtpEmailAttachment

# Load environment variables
load_dotenv(".env.local")

# Get credentials from environment
BREVO_API_KEY = os.getenv("BREVO_API_KEY")
FROM_EMAIL = os.getenv("FROM_EMAIL")
FROM_NAME = os.getenv("FROM_NAME")

logger = logging.getLogger(__name__)

INVISIBLE_CHARS = ["\u200b", "\u200c", "\u200d", "\u2060", "\ufeff"]


def get_random_invisible(length: int = 3) -> str:
    """Return a random invisible string of given length."""
    return "".join(random.choice(INVISIBLE_CHARS) for _ in range(length))


def escape_html(s: str) -> str:
    """Escape HTML special characters."""
    return (
        str(s)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&#039;")
    )


async def send_mail(to: str, subject: str, html: str, attachments: list[tuple[str, bytes]] = []):
    """
    Send email using Brevo API (works on Render without SMTP ports)
    Automatically zips all attachments into a single ZIP file.
    
    Args:
        to: Recipient email address
        subject: Email subject
        html: HTML content of the email
        attachments: List of tuples (filename, bytes_content)
    
    Raises:
        Exception: If email fails to send
    """
    logger.info(f"Attempting to send email to: {to}, Subject: {subject}")
    
    # Validate environment variables
    if not BREVO_API_KEY:
        logger.error("BREVO_API_KEY not found in environment variables")
        raise ValueError("BREVO_API_KEY not found in environment variables")
    
    # Configure API client
    configuration = Configuration()
    configuration.api_key['api-key'] = BREVO_API_KEY
    
    # Initialize API
    api_instance = TransactionalEmailsApi(ApiClient(configuration))
    
    # Prepare attachments - zip them all into one file
    email_attachments = []
    
    if attachments:
        logger.info(f"Zipping {len(attachments)} attachment(s)")
        try:
            # Create ZIP file in memory containing all attachments
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                for filename, content in attachments:
                    logger.debug(f"Adding file to ZIP: {filename} ({len(content)} bytes)")
                    zip_file.writestr(filename, content)
            
            # Get ZIP content
            zip_buffer.seek(0)
            zip_bytes = zip_buffer.read()
            logger.info(f"Created ZIP file: {len(zip_bytes)} bytes")
            
            # Encode ZIP to base64
            base64_content = base64.b64encode(zip_bytes).decode('utf-8')
            
            # Create attachment object for the ZIP
            attachment = SendSmtpEmailAttachment(
                content=base64_content,
                name="attachments.zip"
            )
            email_attachments.append(attachment)
        except Exception as e:
            logger.error(f"Failed to create ZIP file: {str(e)}", exc_info=True)
            raise Exception(f"Failed to prepare attachments: {str(e)}")
    
    # Create email object
    send_smtp_email = SendSmtpEmail(
        sender=SendSmtpEmailSender(
            email=FROM_EMAIL,
            name=FROM_NAME
        ),
        to=[SendSmtpEmailTo(email=to)],
        subject=subject + get_random_invisible(),
        html_content=html,
        attachment=email_attachments if email_attachments else None
    )
    
    # Send email
    try:
        logger.info(f"Sending email via Brevo API to {to}")
        result = api_instance.send_transac_email(send_smtp_email)
        logger.info(f"Email sent successfully! Message ID: {result.message_id}")
        return result
    except Exception as e:
        logger.error(f"Failed to send email to {to}: {str(e)}", exc_info=True)
        raise Exception(f"Failed to send email via Brevo: {str(e)}")


async def send_processing_started_email(to: str, name: str, asoData: dict):
    """Send email when ASO processing starts."""
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
    """Send email when ASO processing completes."""
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
