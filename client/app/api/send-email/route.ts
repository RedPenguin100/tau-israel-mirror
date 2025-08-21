import { NextResponse } from 'next/server';
import nodemailer from 'nodemailer';
export const runtime = 'nodejs';

// .env
// GMAIL_USER=oncoligo.igem@gmail.com
// GMAIL_APP_PASSWORD=xxxx xxxx xxxx xxxx
// EMAIL_FROM="Oncoligo ASO <oncoligo.igem@gmail.com>"

const gmailUser = process.env.GMAIL_USER!;
const gmailPass = process.env.GMAIL_APP_PASSWORD || '';
const fromHeader = process.env.EMAIL_FROM || `Oncoligo ASO <${gmailUser}>`;

const transporter = nodemailer.createTransport({
  host: 'smtp.gmail.com',
  port: 465,           // or 587 with secure: false
  secure: true,        // true for 465, false for 587
  auth: {
    user: gmailUser,   // MUST be the mailbox, not the pretty header
    pass: gmailPass,
  },
});

export async function POST(request: Request) {
  try {
    const { email, name, type, asoData } = await request.json();

    if (!email || !name || !type)
      return NextResponse.json({ error: 'Missing required fields' }, { status: 400 });

    if (!gmailPass)
      return NextResponse.json({ error: 'GMAIL_APP_PASSWORD not set on server' }, { status: 500 });

    if (type === 'processing_started') {
      await sendProcessingStartedEmail(email, name);
    } else if (type === 'processing_completed') {
      await sendProcessingCompletedEmail(email, name, asoData);
    } else {
      return NextResponse.json({ error: 'Unknown email type' }, { status: 400 });
    }

    return NextResponse.json({ success: true });
  } catch (error) {
    console.error('Error sending email:', error);
    return NextResponse.json({ error: 'Failed to send email' }, { status: 500 });
  }
}

async function sendProcessingStartedEmail(to: string, name: string) {
  const subject = 'Your ASO analysis is being processed';
  const html = `
    <div style="font-family: Arial, sans-serif; line-height: 1.6;">
      <h2>Hello ${escapeHtml(name)},</h2>
      <p>We've received your ASO analysis request and processing has begun.</p>
      <p>You'll receive another email with your results once the analysis is complete.</p>
      <p style="color:#555; font-size: 12px;">This message was sent from Oncoligo ASO Generator.</p>
    </div>
  `;
  await transporter.sendMail({ from: fromHeader, to, subject, html });
}

async function sendProcessingCompletedEmail(to: string, name: string, asoData: any = {}) {
  const subject = 'Your ASO analysis is complete!';
  const { organismFile, geneSequence, numericParams } = asoData || {};
  const genePreview =
    typeof geneSequence === 'string'
      ? `${geneSequence.slice(0, 60)}${geneSequence.length > 60 ? '…' : ''}`
      : 'N/A';
  const volume = numericParams?.ASO_volume ?? 'N/A';
  const period = numericParams?.period_of_treatment ?? 'N/A';

  const html = `
    <div style="font-family: Arial, sans-serif; line-height: 1.6;">
      <h2>Hello ${escapeHtml(name)},</h2>
      <p>Your ASO analysis has been completed successfully!</p>
      <h3>Summary</h3>
      <ul>
        <li><strong>Organism file</strong>: ${escapeHtml(String(organismFile || 'N/A'))}</li>
        <li><strong>Gene sequence preview</strong>: <code>${escapeHtml(genePreview)}</code></li>
        <li><strong>ASO volume</strong>: ${escapeHtml(String(volume))}</li>
        <li><strong>Treatment period</strong>: ${escapeHtml(String(period))} days</li>
      </ul>
      <p style="color:#555; font-size: 12px;">This message was sent from Oncoligo ASO Generator.</p>
    </div>
  `;
  await transporter.sendMail({ from: fromHeader, to, subject, html });
}

function escapeHtml(input: any) {
  return String(input)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#039;');
}
