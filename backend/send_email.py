import os
import base64
import io
import zipfile
import logging
import random
import json
from dotenv import load_dotenv
from brevo_python import Configuration, ApiClient, TransactionalEmailsApi
from brevo_python import SendSmtpEmail, SendSmtpEmailSender, SendSmtpEmailTo, SendSmtpEmailAttachment
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import pandas as pd
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, A4
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer, PageBreak
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib.enums import TA_CENTER, TA_LEFT

# Load environment variables
load_dotenv(".env.local")

# Get credentials from environment
BREVO_API_KEY = os.getenv("BREVO_API_KEY")
FROM_EMAIL = os.getenv("FROM_EMAIL")
FROM_NAME = os.getenv("FROM_NAME")

logger = logging.getLogger(__name__)

INVISIBLE_CHARS = ["\u200b", "\u200c", "\u200d", "\u2060", "\ufeff"]

# Modification pattern explanation
MODIFICATION_EXPLANATION = """
Modification Pattern Guide:
• L = LNA (Locked Nucleic Acid) - Enhanced binding affinity
• M = MOE (2'-O-Methoxyethyl) - Improved nuclease resistance
• d = DNA (Unmodified deoxyribonucleotide)
• All ASO sequences have phosphorothioate (PS) backbone for stability

Example: LLLddddddddddLLL = LNA-DNA-LNA gapmer with PS backbone
Example: MMMMMddddddddddMMMMM = MOE-DNA-MOE gapmer with PS backbone
"""


def create_sbol_format(df, full_mrna_seq, gene_name):
    """
    Create SBOL (Synthetic Biology Open Language) format output.
    Simplified SBOL v2 representation for ASO sequences.
    """
    # Handle None or empty gene names
    if not gene_name or gene_name == "None":
        gene_name = "Selected_Gene"
    
    sbol_lines = []
    sbol_lines.append("<?xml version='1.0' encoding='UTF-8'?>")
    sbol_lines.append("<rdf:RDF xmlns:rdf='http://www.w3.org/1999/02/22-rdf-syntax-ns#'")
    sbol_lines.append("         xmlns:sbol='http://sbols.org/v2#'>")
    sbol_lines.append("  <!-- Modification Pattern Guide: L=LNA, M=MOE, d=DNA (unmodified), All sequences have PS backbone -->")
    
    for idx, row in df.iterrows():
        seq = row['Sequence']
        aso_id = f"ASO_{idx}"
        
        sbol_lines.append(f"  <sbol:ComponentDefinition rdf:about='#{aso_id}'>")
        sbol_lines.append(f"    <sbol:displayId>{aso_id}</sbol:displayId>")
        sbol_lines.append(f"    <sbol:name>{gene_name}_ASO_{idx}</sbol:name>")
        sbol_lines.append("    <sbol:type rdf:resource='http://www.biopax.org/release/biopax-level3.owl#DnaRegion'/>")
        sbol_lines.append("    <sbol:role rdf:resource='http://identifiers.org/so/SO:0000336'/>") # Antisense RNA
        sbol_lines.append(f"    <sbol:sequence rdf:resource='#seq_{aso_id}'/>")
        
        if 'mod_pattern' in row:
            sbol_lines.append(f"    <sbol:description>Modification: {row['mod_pattern']} (L=LNA, M=MOE, d=DNA, PS backbone)</sbol:description>")
        
        sbol_lines.append("  </sbol:ComponentDefinition>")
        
        sbol_lines.append(f"  <sbol:Sequence rdf:about='#seq_{aso_id}'>")
        sbol_lines.append(f"    <sbol:elements>{seq}</sbol:elements>")
        sbol_lines.append("    <sbol:encoding rdf:resource='http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html'/>")
        sbol_lines.append("  </sbol:Sequence>")
    
    sbol_lines.append("</rdf:RDF>")
    return "\n".join(sbol_lines)


def create_genbank_format(df, full_mrna_seq, gene_name, organism_name="Unknown"):
    """
    Create GenBank format file with ASO annotations on the mRNA sequence.
    """
    # Handle None or empty gene names
    if not gene_name or gene_name == "None":
        gene_name = "Selected_Gene"
    
    gb_records = []
    
    for idx, row in df.iterrows():
        seq = row['Sequence']
        sense_start = int(row.get('sense_start', 0)) if 'sense_start' in row else 0
        
        record = SeqRecord(
            Seq(seq),
            id=f"ASO_{idx}",
            name=f"{gene_name}_ASO",
            description=f"Antisense Oligonucleotide for {gene_name} at position {sense_start}"
        )
        
        record.annotations["molecule_type"] = "DNA"
        record.annotations["topology"] = "linear"
        record.annotations["organism"] = organism_name
        record.annotations["note"] = "L=LNA, M=MOE, d=DNA (unmodified), All sequences have phosphorothioate (PS) backbone"
        
        if 'mod_pattern' in row:
            record.annotations["modification"] = row['mod_pattern']
        if 'gc_content' in row:
            record.annotations["gc_content"] = f"{row['gc_content']:.2f}"
        
        gb_records.append(record)
    
    gb_io = StringIO()
    SeqIO.write(gb_records, gb_io, "genbank")
    return gb_io.getvalue()


def create_visualization_html(df, full_mrna_seq, gene_name):
    """
    Create an interactive HTML visualization showing ASO positions on mRNA.
    """
    # Handle None or empty gene names
    if not gene_name or gene_name == "None":
        gene_name = "Selected_Gene"
    
    html_parts = []
    html_parts.append("""
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>ASO Visualization</title>
    <style>
        body { font-family: 'Segoe UI', Arial, sans-serif; margin: 20px; background: #f5f5f5; }
        .container { max-width: 1200px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
        h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
        .info-box { background: #e8f4f8; border-left: 4px solid #3498db; padding: 15px; margin: 20px 0; border-radius: 4px; }
        .info-box h3 { margin-top: 0; color: #2c3e50; }
        .info-box code { background: white; padding: 2px 6px; border-radius: 3px; font-family: 'Courier New', monospace; }
        .website-link { background: #3498db; color: white; padding: 12px 24px; text-decoration: none; border-radius: 5px; display: inline-block; margin: 15px 0; font-weight: bold; transition: background 0.3s; }
        .website-link:hover { background: #2980b9; }
        .mrna-viz { margin: 30px 0; position: relative; }
        .mrna-bar { height: 40px; background: linear-gradient(to right, #3498db, #2ecc71); border-radius: 5px; position: relative; margin: 20px 0; }
        .aso-marker { position: absolute; height: 100%; background: rgba(231, 76, 60, 0.7); border: 2px solid #c0392b; border-radius: 3px; cursor: pointer; transition: all 0.3s; }
        .aso-marker:hover { background: rgba(231, 76, 60, 0.9); transform: translateY(-2px); box-shadow: 0 4px 8px rgba(0,0,0,0.2); }
        .aso-label { position: absolute; top: -25px; left: 50%; transform: translateX(-50%); font-size: 11px; font-weight: bold; color: #c0392b; white-space: nowrap; }
        .legend { margin: 20px 0; padding: 15px; background: #ecf0f1; border-radius: 5px; }
        .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }
        .stat-card { background: #3498db; color: white; padding: 15px; border-radius: 5px; text-align: center; }
        .stat-value { font-size: 24px; font-weight: bold; }
        .stat-label { font-size: 12px; opacity: 0.9; }
        table { width: 100%; border-collapse: collapse; margin: 20px 0; }
        th, td { padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }
        th { background: #34495e; color: white; font-weight: 600; }
        tr:hover { background: #f8f9fa; }
        .sequence { word-break: break-all; background: #ecf0f1; padding: 10px; border-radius: 3px; }
        .tooltip { position: absolute; background: #2c3e50; color: white; padding: 10px; border-radius: 5px; font-size: 12px; pointer-events: none; z-index: 1000; display: none; box-shadow: 0 4px 12px rgba(0,0,0,0.3); }
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 ASO Analysis Results: """ + gene_name + """</h1>
        
        <div class="info-box">
            <h3>🔬 Modification Pattern Guide</h3>
            <p><strong>L</strong> = LNA (Locked Nucleic Acid) - Enhanced binding affinity<br>
            <strong>M</strong> = MOE (2'-O-Methoxyethyl) - Improved nuclease resistance<br>
            <strong>d</strong> = DNA (Unmodified deoxyribonucleotide)<br>
            <strong>All ASO sequences have phosphorothioate (PS) backbone for stability</strong></p>
            <p style="margin-top: 10px;">
            <strong>Example:</strong> <code>LLLddddddddddLLL</code> = LNA-DNA-LNA gapmer with PS backbone<br>
            <strong>Example:</strong> <code>MMMMMddddddddddMMMMM</code> = MOE-DNA-MOE gapmer with PS backbone
            </p>
        </div>
        
        <a href="https://2025.igem.wiki/tau-israel/model" target="_blank" class="website-link">
            📚 Learn More About Our ASO Design Model
        </a>
        
        <div class="stats-grid">
            <div class="stat-card">
                <div class="stat-value">""" + str(len(df)) + """</div>
                <div class="stat-label">Total ASO Sequences</div>
            </div>
            <div class="stat-card" style="background: #2ecc71;">
                <div class="stat-value">""" + str(len(full_mrna_seq)) + """</div>
                <div class="stat-label">mRNA Length (nt)</div>
            </div>
""")
    
    if 'gc_content' in df.columns:
        avg_gc = df['gc_content'].mean()
        html_parts.append(f"""
            <div class="stat-card" style="background: #e74c3c;">
                <div class="stat-value">{avg_gc:.1f}%</div>
                <div class="stat-label">Avg GC Content</div>
            </div>
""")
    
    html_parts.append("""
        </div>
        
        <h2>ASO Binding Sites on mRNA</h2>
        <div class="mrna-viz">
            <div style="display: flex; justify-content: space-between; margin-bottom: 5px; font-size: 12px; color: #7f8c8d;">
                <span>5' Start</span>
                <span>3' End (""" + str(len(full_mrna_seq)) + """ nt)</span>
            </div>
            <div class="mrna-bar" id="mrna-bar">
""")
    
    mrna_len = len(full_mrna_seq)
    for idx, row in df.iterrows():
        if 'sense_start' in row and 'Sequence' in row:
            start = int(row['sense_start'])
            seq_len = len(row['Sequence'])
            left_pct = (start / mrna_len) * 100
            width_pct = (seq_len / mrna_len) * 100
            
            tooltip_info = f"ASO {idx}"
            if 'mod_pattern' in row:
                tooltip_info += f"<br>Mod: {row['mod_pattern']}"
            if 'gc_content' in row:
                tooltip_info += f"<br>GC: {row['gc_content']:.1f}%"
            
            html_parts.append(f"""
                <div class="aso-marker" style="left: {left_pct}%; width: {width_pct}%;" 
                     onmouseover="showTooltip(event, '{tooltip_info}')" 
                     onmouseout="hideTooltip()">
                    <div class="aso-label">ASO {idx}</div>
                </div>
""")
    
    html_parts.append("""
            </div>
        </div>
        
        <div class="legend">
            <strong>Legend:</strong> 
            <span style="display: inline-block; width: 20px; height: 20px; background: linear-gradient(to right, #3498db, #2ecc71); border-radius: 3px; vertical-align: middle; margin: 0 5px;"></span> mRNA Sequence
            <span style="display: inline-block; width: 20px; height: 20px; background: rgba(231, 76, 60, 0.7); border: 2px solid #c0392b; border-radius: 3px; vertical-align: middle; margin: 0 5px 0 15px;"></span> ASO Binding Site
        </div>
        
        <h2>Detailed ASO Information</h2>
        <table>
            <thead>
                <tr>
                    <th>ID</th>
                    <th>Position</th>
                    <th>Length</th>
""")
    
    if 'mod_pattern' in df.columns:
        html_parts.append("<th>Modification</th>")
    if 'gc_content' in df.columns:
        html_parts.append("<th>GC%</th>")
    if 'exp_ps_hybr' in df.columns:
        html_parts.append("<th>Exp PS Hybr</th>")
    
    html_parts.append("""
                    <th>Sequence</th>
                </tr>
            </thead>
            <tbody>
""")
    
    for idx, row in df.iterrows():
        start = int(row.get('sense_start', 0)) if 'sense_start' in row else 'N/A'
        seq_len = len(row['Sequence'])
        
        html_parts.append(f"""
                <tr>
                    <td><strong>ASO {idx}</strong></td>
                    <td>{start}</td>
                    <td>{seq_len} nt</td>
""")
        
        if 'mod_pattern' in df.columns:
            html_parts.append(f"<td>{row.get('mod_pattern', 'N/A')}</td>")
        if 'gc_content' in df.columns:
            html_parts.append(f"<td>{row.get('gc_content', 0):.1f}%</td>")
        if 'exp_ps_hybr' in df.columns:
            html_parts.append(f"<td>{row.get('exp_ps_hybr', 'N/A'):.2f}</td>")
        
        html_parts.append(f"""
                    <td class="sequence">{row['Sequence']}</td>
                </tr>
""")
    
    html_parts.append("""
            </tbody>
        </table>
        
        <div style="text-align: center; margin-top: 40px; padding: 20px; background: #f8f9fa; border-radius: 5px;">
            <p style="margin-bottom: 15px;">For more information about our ASO design methodology:</p>
            <a href="https://2025.igem.wiki/tau-israel/model" target="_blank" class="website-link">
                🌐 Visit Our Model Documentation
            </a>
        </div>
    </div>
    
    <div class="tooltip" id="tooltip"></div>
    
    <script>
        function showTooltip(event, text) {
            const tooltip = document.getElementById('tooltip');
            tooltip.innerHTML = text;
            tooltip.style.display = 'block';
            tooltip.style.left = event.pageX + 10 + 'px';
            tooltip.style.top = event.pageY - 30 + 'px';
        }
        
        function hideTooltip() {
            document.getElementById('tooltip').style.display = 'none';
        }
    </script>
</body>
</html>
""")
    
    return "".join(html_parts)


def create_pdf_report(df, full_mrna_seq, gene_name, organism_name="Unknown"):
    """
    Create a comprehensive PDF report with ASO analysis results.
    """
    # Handle None or empty gene names
    if not gene_name or gene_name == "None":
        gene_name = "Selected Gene"
    
    buffer = io.BytesIO()
    # Reduce margins to make more room for the table
    doc = SimpleDocTemplate(buffer, pagesize=letter, 
                           topMargin=0.4*inch, bottomMargin=0.4*inch,
                           leftMargin=0.4*inch, rightMargin=0.4*inch)
    story = []
    styles = getSampleStyleSheet()
    
    # Custom styles
    title_style = ParagraphStyle(
        'CustomTitle',
        parent=styles['Heading1'],
        fontSize=24,
        textColor=colors.HexColor('#2c3e50'),
        spaceAfter=30,
        alignment=TA_CENTER
    )
    
    heading_style = ParagraphStyle(
        'CustomHeading',
        parent=styles['Heading2'],
        fontSize=16,
        textColor=colors.HexColor('#3498db'),
        spaceAfter=12,
        spaceBefore=12
    )
    
    info_style = ParagraphStyle(
        'InfoBox',
        parent=styles['BodyText'],
        fontSize=9,
        textColor=colors.HexColor('#2c3e50'),
        leftIndent=10,
        rightIndent=10,
        spaceAfter=6
    )
    
    # Title
    story.append(Paragraph(f"ASO Analysis Report: {gene_name}", title_style))
    story.append(Spacer(1, 0.1*inch))
    
    # Modification guide box
    story.append(Paragraph("Modification Pattern Guide", heading_style))
    mod_text = """
    <b>L</b> = LNA (Locked Nucleic Acid) - Enhanced binding affinity<br/>
    <b>M</b> = MOE (2'-O-Methoxyethyl) - Improved nuclease resistance<br/>
    <b>d</b> = DNA (Unmodified deoxyribonucleotide)<br/>
    <b>All ASO sequences have phosphorothioate (PS) backbone for stability</b><br/>
    <br/>
    <b>Example:</b> LLLddddddddddLLL = LNA-DNA-LNA gapmer with PS backbone<br/>
    <b>Example:</b> MMMMMddddddddddMMMMM = MOE-DNA-MOE gapmer with PS backbone
    """
    story.append(Paragraph(mod_text, info_style))
    story.append(Spacer(1, 0.15*inch))
    
    # Website link
    website_text = """
    For more information about our ASO design methodology, visit:<br/>
    <link href="https://2025.igem.wiki/tau-israel/model" color="blue">
    https://2025.igem.wiki/tau-israel/model
    </link>
    """
    story.append(Paragraph(website_text, info_style))
    story.append(Spacer(1, 0.2*inch))
    
    # Summary section
    story.append(Paragraph("Summary", heading_style))
    summary_data = [
        ["Parameter", "Value"],
        ["Gene Name", gene_name],
        ["Organism", organism_name],
        ["Total ASO Sequences", str(len(df))],
        ["mRNA Length", f"{len(full_mrna_seq)} nucleotides"],
    ]
    
    if 'gc_content' in df.columns:
        avg_gc = df['gc_content'].mean()
        summary_data.append(["Average GC Content", f"{avg_gc:.2f}%"])
    
    summary_table = Table(summary_data, colWidths=[2.5*inch, 4.5*inch])
    summary_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#34495e')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 12),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
    ]))
    story.append(summary_table)
    story.append(Spacer(1, 0.3*inch))
    
    # ASO Details Table
    story.append(Paragraph("ASO Sequences", heading_style))
    
    # Prepare table data
    table_data = [["ID", "Position", "Length", "Sequence"]]
    
    # Add optional columns
    optional_cols = []
    if 'mod_pattern' in df.columns:
        table_data[0].insert(-1, "Modification")
        optional_cols.append('mod_pattern')
    if 'gc_content' in df.columns:
        table_data[0].insert(-1, "GC%")
        optional_cols.append('gc_content')
    
    for idx, row in df.iterrows():
        row_data = [
            f"ASO {idx}",
            str(int(row.get('sense_start', 0))) if 'sense_start' in row else 'N/A',
            str(len(row['Sequence'])),
        ]
        
        for col in optional_cols:
            if col == 'gc_content':
                row_data.append(f"{row[col]:.1f}%")
            else:
                row_data.append(str(row[col]))
        
        # Format sequence - don't truncate, use full sequence (up to 20nt fits well)
        seq = row['Sequence']
        if len(seq) > 20:
            seq = seq[:20] + "..."
        row_data.append(seq)
        
        table_data.append(row_data)
    
    # Calculate column widths dynamically - WIDER ID COLUMN (0.75" for ASO1xx = 10 chars)
    # Available width: 8.5" - 0.8" (margins) = 7.7"
    num_cols = len(table_data[0])
    if num_cols == 4:  # ID, Position, Length, Sequence
        col_widths = [0.75*inch, 0.9*inch, 0.7*inch, 5.35*inch]
    elif num_cols == 5:  # + Modification or GC%
        col_widths = [0.75*inch, 0.8*inch, 0.6*inch, 2.6*inch, 2.95*inch]
    elif num_cols == 6:  # + Both Modification and GC%
        col_widths = [0.75*inch, 0.7*inch, 0.5*inch, 2.4*inch, 0.75*inch, 2.6*inch]
    else:  # More columns
        # Distribute space: prioritize sequence and modification columns
        fixed_width = 1.95*inch  # ID, Position, Length (ID is now 0.75")
        seq_width = 2.4*inch  # Sequence column
        mod_width = 2.4*inch  # Modification column (if present)
        remaining = 7.7*inch - fixed_width - seq_width
        var_cols = num_cols - 4  # subtract fixed cols and sequence
        
        if 'mod_pattern' in df.columns:
            remaining -= mod_width
            var_cols -= 1
        
        var_width = (remaining / var_cols) if var_cols > 0 else 0.8*inch
        
        # Build column widths list
        col_widths = [0.75*inch, 0.7*inch, 0.5*inch]
        for i in range(3, num_cols - 1):
            if i == 3 and 'mod_pattern' in df.columns:
                col_widths.append(mod_width)
            else:
                col_widths.append(var_width)
        col_widths.append(seq_width)  # Last column is always sequence
    
    aso_table = Table(table_data, colWidths=col_widths)
    aso_table.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#3498db')),
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
        ('ALIGN', (0, 0), (-1, -1), 'LEFT'),
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, 0), 10),
        ('FONTSIZE', (0, 1), (-1, -1), 8),
        ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
        ('GRID', (0, 0), (-1, -1), 1, colors.black),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.lightgrey]),
    ]))
    story.append(aso_table)
    
    # Build PDF
    doc.build(story)
    buffer.seek(0)
    return buffer.read()


def prepare_data_before_sending(result_from_foo):
    """
    Prepare comprehensive data files for biologists and bioinformaticians.
    
    Args:
        result_from_foo: Tuple of (DataFrame, full_mRNA_sequence, gene_name, organism_name)
    
    Returns:
        List of tuples (filename, bytes_content) for all generated files
    """
    # Unpack result
    if len(result_from_foo) == 4:
        df, full_mrna_seq, gene_name, organism_name = result_from_foo
    elif len(result_from_foo) == 3:
        df, full_mrna_seq, gene_name = result_from_foo
        organism_name = "Unknown"
    else:
        df, full_mrna_seq = result_from_foo
        gene_name = "Selected_Gene"
        organism_name = "Unknown"
    
    # Handle None or empty gene names
    if not gene_name or gene_name == "None":
        gene_name = "Selected_Gene"
    
    files = []
    
    # 1. FASTA format (standard molecular biology format)
    try:
        fasta_records = []
        fasta_header = f"# Modification Guide: L=LNA, M=MOE, d=DNA (unmodified), All sequences have PS backbone\n"
        for idx, row in df.iterrows():
            description_parts = [f"ASO_{idx}"]
            if 'mod_pattern' in df.columns:
                description_parts.append(f"mod={row['mod_pattern']}")
            if 'sense_start' in df.columns:
                description_parts.append(f"pos={int(row['sense_start'])}")
            if 'gc_content' in df.columns:
                description_parts.append(f"GC={row['gc_content']:.1f}%")
            
            record = SeqRecord(
                Seq(row['Sequence']),
                id=f"ASO_{idx}",
                description=" | ".join(description_parts)
            )
            fasta_records.append(record)
        
        fasta_io = StringIO()
        fasta_io.write(fasta_header)
        SeqIO.write(fasta_records, fasta_io, "fasta")
        files.append(("aso_sequences.fasta", fasta_io.getvalue().encode()))
    except Exception as e:
        logger.error(f"Error creating FASTA: {e}")
    
    # 2. CSV format (for easy data analysis)
    try:
        csv_io = StringIO()
        csv_io.write("# Modification Guide: L=LNA, M=MOE, d=DNA (unmodified), All sequences have PS backbone\n")
        df.to_csv(csv_io, index=True, index_label="ASO_ID")
        files.append(("aso_data.csv", csv_io.getvalue().encode()))
    except Exception as e:
        logger.error(f"Error creating CSV: {e}")
    
    # 3. JSON format (for programmatic access)
    try:
        json_data = {
            "modification_guide": {
                "L": "LNA (Locked Nucleic Acid) - Enhanced binding affinity",
                "M": "MOE (2'-O-Methoxyethyl) - Improved nuclease resistance",
                "d": "DNA (Unmodified deoxyribonucleotide)",
                "backbone": "All sequences have phosphorothioate (PS) backbone"
            },
            "gene_name": gene_name,
            "organism": organism_name,
            "mrna_length": len(full_mrna_seq),
            "mrna_sequence": full_mrna_seq,
            "aso_count": len(df),
            "aso_sequences": df.to_dict(orient='records'),
            "reference": "https://2025.igem.wiki/tau-israel/model"
        }
        json_str = json.dumps(json_data, indent=2)
        files.append(("aso_data.json", json_str.encode()))
    except Exception as e:
        logger.error(f"Error creating JSON: {e}")
    
    # 4. SBOL format (Synthetic Biology standard)
    try:
        sbol_content = create_sbol_format(df, full_mrna_seq, gene_name)
        files.append(("aso_sequences.sbol", sbol_content.encode()))
    except Exception as e:
        logger.error(f"Error creating SBOL: {e}")
    
    # 5. GenBank format
    try:
        genbank_content = create_genbank_format(df, full_mrna_seq, gene_name, organism_name)
        files.append(("aso_sequences.gb", genbank_content.encode()))
    except Exception as e:
        logger.error(f"Error creating GenBank: {e}")
    
    # 6. HTML Visualization
    try:
        html_content = create_visualization_html(df, full_mrna_seq, gene_name)
        files.append(("aso_visualization.html", html_content.encode()))
    except Exception as e:
        logger.error(f"Error creating HTML visualization: {e}")
    
    # 7. PDF Report
    try:
        pdf_content = create_pdf_report(df, full_mrna_seq, gene_name, organism_name)
        files.append(("aso_report.pdf", pdf_content))
    except Exception as e:
        logger.error(f"Error creating PDF report: {e}")
    
    # 8. Full mRNA sequence as separate FASTA (only if we have a real gene name)
    if gene_name and gene_name != "None" and gene_name != "Selected_Gene":
        try:
            mrna_record = SeqRecord(
                Seq(full_mrna_seq),
                id=gene_name,
                description=f"Full mRNA sequence for {gene_name} ({organism_name})"
            )
            mrna_io = StringIO()
            SeqIO.write([mrna_record], mrna_io, "fasta")
            files.append(("full_mrna.fasta", mrna_io.getvalue().encode()))
        except Exception as e:
            logger.error(f"Error creating mRNA FASTA: {e}")
    
    # 9. Statistics summary text file
    try:
        stats_lines = []
        stats_lines.append(f"ASO Analysis Summary for {gene_name}")
        stats_lines.append("=" * 60)
        stats_lines.append("")
        stats_lines.append("MODIFICATION PATTERN GUIDE:")
        stats_lines.append("L = LNA (Locked Nucleic Acid) - Enhanced binding affinity")
        stats_lines.append("M = MOE (2'-O-Methoxyethyl) - Improved nuclease resistance")
        stats_lines.append("d = DNA (Unmodified deoxyribonucleotide)")
        stats_lines.append("All ASO sequences have phosphorothioate (PS) backbone")
        stats_lines.append("")
        stats_lines.append("Example: LLLddddddddddLLL = LNA-DNA-LNA gapmer with PS backbone")
        stats_lines.append("Example: MMMMMddddddddddMMMMM = MOE-DNA-MOE gapmer with PS backbone")
        stats_lines.append("")
        stats_lines.append("=" * 60)
        stats_lines.append(f"Organism: {organism_name}")
        stats_lines.append(f"Gene: {gene_name}")
        stats_lines.append(f"mRNA Length: {len(full_mrna_seq)} nucleotides")
        stats_lines.append(f"Total ASO Sequences: {len(df)}")
        stats_lines.append("")
        
        if 'gc_content' in df.columns:
            stats_lines.append(f"Average GC Content: {df['gc_content'].mean():.2f}%")
            stats_lines.append(f"GC Content Range: {df['gc_content'].min():.2f}% - {df['gc_content'].max():.2f}%")
        
        if 'sense_start' in df.columns:
            stats_lines.append(f"Position Range: {df['sense_start'].min():.0f} - {df['sense_start'].max():.0f}")
        
        if 'exp_ps_hybr' in df.columns:
            stats_lines.append(f"Average Exp PS Hybr: {df['exp_ps_hybr'].mean():.3f}")
        
        stats_lines.append("")
        stats_lines.append("For more information about our ASO design methodology:")
        stats_lines.append("https://2025.igem.wiki/tau-israel/model")
        
        stats_content = "\n".join(stats_lines)
        files.append(("summary_statistics.txt", stats_content.encode()))
    except Exception as e:
        logger.error(f"Error creating statistics summary: {e}")
    
    logger.info(f"Prepared {len(files)} files for email attachment")
    return files


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
                name="aso_analysis_results.zip"
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


async def send_processing_completed_email(to: str, name: str, asoData: dict, files: list[tuple[str, bytes]] | None):
    """Send email when ASO processing completes with comprehensive attachments."""
    subject = "Your ASO analysis is complete!"

    organism_file = asoData.get("organismFile", "N/A")
    gene_name = asoData.get("geneName", "N/A")
    gene_seq = asoData.get("geneSequence") or ""
    gene_preview = gene_seq[:60] + ("…" if len(gene_seq) > 60 else "")
    top_k = asoData.get("top_k", "N/A")
    include_breakdown = "Yes" if asoData.get("includeFeatureBreakdown") else "No"

    # Determine which files were generated
    file_list_html = ""
    if files:
        file_list_html = "<h3>📦 Attached Files</h3><ul>"
        for filename, _ in files:
            file_desc = {
                "aso_sequences.fasta": "FASTA format - Standard sequence file",
                "aso_data.csv": "CSV format - Spreadsheet-compatible data",
                "aso_data.json": "JSON format - Programmatic access",
                "aso_sequences.sbol": "SBOL format - Synthetic biology standard",
                "aso_sequences.gb": "GenBank format - Annotated sequences",
                "aso_visualization.html": "HTML Visualization - Interactive viewer",
                "aso_report.pdf": "PDF Report - Comprehensive analysis",
                "full_mrna.fasta": "mRNA sequence - Full target sequence",
                "summary_statistics.txt": "Statistics - Quick summary"
            }
            desc = file_desc.get(filename, "Additional data file")
            file_list_html += f"<li><strong>{filename}</strong> - {desc}</li>"
        file_list_html += "</ul>"

    html = f"""
    <div style="font-family: Arial, sans-serif; line-height: 1.6;">
      <h2>Hello {escape_html(name)},</h2>
      <p>Your ASO analysis has been completed successfully! 🎉</p>
      
      <h3>Summary</h3>
      <ul>
        <li><strong>Organism file</strong>: {escape_html(str(organism_file))}</li>
        <li><strong>Gene name</strong>: {escape_html(str(gene_name))}</li>
        <li><strong>Gene sequence preview</strong>: <code>{escape_html(gene_preview or "N/A")}</code></li>
        <li><strong>Top-k results</strong>: {escape_html(str(top_k))}</li>
        <li><strong>Detailed analysis</strong>: {escape_html(include_breakdown)}</li>
      </ul>
      
      <div style="background: #e8f4f8; border-left: 4px solid #3498db; padding: 15px; margin: 20px 0;">
        <h3 style="margin-top: 0;">🔬 Modification Pattern Guide</h3>
        <p><strong>L</strong> = LNA (Locked Nucleic Acid) - Enhanced binding affinity<br>
        <strong>M</strong> = MOE (2'-O-Methoxyethyl) - Improved nuclease resistance<br>
        <strong>d</strong> = DNA (Unmodified deoxyribonucleotide)<br>
        <strong>All ASO sequences have phosphorothioate (PS) backbone for stability</strong></p>
        <p><strong>Example:</strong> <code>LLLddddddddddLLL</code> = LNA-DNA-LNA gapmer with PS backbone<br>
        <strong>Example:</strong> <code>MMMMMddddddddddMMMMM</code> = MOE-DNA-MOE gapmer with PS backbone</p>
      </div>
      
      {file_list_html}
      
      <p>All files are included in the attached <strong>aso_analysis_results.zip</strong> file.</p>
      
      <h3>📖 Quick Start Guide</h3>
      <ul>
        <li><strong>For visualization</strong>: Open <code>aso_visualization.html</code> in your web browser</li>
        <li><strong>For detailed report</strong>: Open <code>aso_report.pdf</code></li>
        <li><strong>For data analysis</strong>: Use <code>aso_data.csv</code> in Excel/Python/R</li>
        <li><strong>For molecular biology tools</strong>: Use <code>aso_sequences.fasta</code> or <code>aso_sequences.gb</code></li>
      </ul>
      
      <div style="text-align: center; margin: 30px 0; padding: 20px; background: #f8f9fa; border-radius: 5px;">
        <p style="margin-bottom: 15px;"><strong>For more information about our ASO design methodology:</strong></p>
        <a href="https://2025.igem.wiki/tau-israel/model" 
           style="background: #3498db; color: white; padding: 12px 24px; text-decoration: none; 
                  border-radius: 5px; display: inline-block; font-weight: bold;">
          🌐 Visit Our Model Documentation
        </a>
      </div>
      
      <p style="margin-top: 20px;">Thank you for using Oncoligo ASO Generator!</p>
      <p style="color:#555; font-size: 12px;">This message was sent from Oncoligo ASO Generator.</p>
    </div>
    """

    await send_mail(to, subject, html, files or [])