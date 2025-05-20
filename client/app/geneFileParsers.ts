/**
 * Parse sequence files (GenBank, FASTA, JBEI, SnapGene, SBOL) or accession IDs (NCBI, iGEM) to a simple, common format
 *
 * adapted from: https://github.com/Lattice-Automation/seqparse
 **/

/* Seq is a single parsed sequence from a file or accession. */
export interface Seq {
  /** name of the sequence */
  name: string;
  /** the sequence */
  seq: string;
  /** type of sequence. Inferred from the seq's symbols */
  type: "dna" | "rna" | "aa" | "unknown";
}

/**
 * mapping the 64 standard codons to amino acids
 * no synth AA's
 *
 * adapted from: "https://github.com/keithwhor/NtSeq/blob/master/lib/nt.js
 **/
const codon2AA = {
  AAA: "K",
  AAC: "N",
  AAG: "K",
  AAT: "N",
  ACA: "T",
  ACC: "T",
  ACG: "T",
  ACT: "T",
  AGA: "R",
  AGC: "S",
  AGG: "R",
  AGT: "S",
  ATA: "I",
  ATC: "I",
  ATG: "M",
  ATT: "I",
  CAA: "Q",
  CAC: "H",
  CAG: "Q",
  CAT: "H",
  CCA: "P",
  CCC: "P",
  CCG: "P",
  CCT: "P",
  CGA: "R",
  CGC: "R",
  CGG: "R",
  CGT: "R",
  CTA: "L",
  CTC: "L",
  CTG: "L",
  CTT: "L",
  GAA: "E",
  GAC: "D",
  GAG: "E",
  GAT: "D",
  GCA: "A",
  GCC: "A",
  GCG: "A",
  GCT: "A",
  GGA: "G",
  GGC: "G",
  GGG: "G",
  GGT: "G",
  GTA: "V",
  GTC: "V",
  GTG: "V",
  GTT: "V",
  TAA: "*",
  TAC: "Y",
  TAG: "*",
  TAT: "Y",
  TCA: "S",
  TCC: "S",
  TCG: "S",
  TCT: "S",
  TGA: "*",
  TGC: "C",
  TGG: "W",
  TGT: "C",
  TTA: "L",
  TTC: "F",
  TTG: "L",
  TTT: "F",
};

// from http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
const comp = {
  A: "T",
  B: "V",
  C: "G",
  D: "H",
  G: "C",
  H: "D",
  K: "M",
  M: "K",
  N: "N",
  R: "Y",
  S: "S",
  T: "A",
  U: "A",
  V: "B",
  W: "W",
  X: "X",
  Y: "R",
  a: "t",
  b: "v",
  c: "g",
  d: "h",
  g: "c",
  h: "d",
  k: "m",
  m: "k",
  n: "n",
  r: "y",
  s: "s",
  t: "a",
  u: "a",
  v: "b",
  w: "w",
  x: "x",
  y: "r",
};

/** Infer the type of a sequence. This only allows a couple wildcard characters so may be overly strict. */
const aminoAcids = Array.from(new Set(Object.values(codon2AA)).values()).join(
  ""
);
const aminoAcidRegex = new RegExp(`^[${aminoAcids}]+$`, "i");
const guessType = (seq: string): "dna" | "rna" | "aa" | "unknown" => {
  if (/^[atgcn.]+$/i.test(seq)) {
    return "dna";
  } else if (/^[augcn.]+$/i.test(seq)) {
    return "rna";
  } else if (aminoAcidRegex.test(seq)) {
    return "aa";
  }
  return "unknown";
};

/**
 * Return the filtered sequence and its complement if its an empty string, return the same for both.
 */
export const complement = (
  origSeq: string
): { compSeq: string; seq: string } => {
  if (!origSeq) {
    return { compSeq: "", seq: "" };
  }

  // filter out unrecognized basepairs and build up the complement
  let seq = "";
  let compSeq = "";
  for (let i = 0, origLength = origSeq.length; i < origLength; i += 1) {
    if (comp[origSeq[i] as keyof typeof comp]) {
      seq += origSeq[i];
      compSeq += comp[origSeq[i] as keyof typeof comp];
    }
  }
  return { compSeq, seq };
};

export const fastaParser = (text: string, fileName: string): Seq[] => {
  // partFactory returns a negative "circular" prop, we assume they're all linear
  if (text.trim().startsWith(">")) {
    return text
      .split(">") // split up if it's a multi-seq FASTA file
      .map((t) => {
        // this starts at the end of the first line, grabs all other characters,
        // and removes any newlines (leaving only the original sequence)
        // sequence "cleaning" happens in complement (we don't support bps other than
        // the most common right now)
        const seq = t.substr(t.indexOf("\n"), t.length).replace(/\s/g, "");

        // the first line contains the name, though there's lots of variability around
        // the information on this line...
        // >MCHU - Calmodulin - Human, rabbit, bovine, rat, and chicken
        const name = t.substring(0, t.search(/\n|\|/)).replace(/\//g, "");

        return {
          name,
          seq,
          type: guessType(seq),
        };
      })
      .filter((p) => p.name && p.seq);
  }

  if (text.trim().startsWith(";")) {
    // it's an old-school style FASTA that's punctuated with semi-colons
    // ;my|NAME
    // ;my comment
    // actGacgata
    const name = text.substring(0, text.search(/\n|\|/)).replace(/\//g, "");
    const newlineBeforeSeq = text.indexOf("\n", text.lastIndexOf(";"));
    const seq = text.substring(newlineBeforeSeq, text.length);
    return [
      {
        name,
        seq,
        type: guessType(seq),
      },
    ];
  }

  // assume that it's a no name FASTA. Ie it's just a file with dna and no header
  // try and get the name from the fileName
  const lastChar = fileName.lastIndexOf(".") || fileName.length;
  const name = fileName.substring(0, lastChar) || "Untitled";
  const seq = text;
  return [
    {
      name,
      seq,
      type: guessType(seq),
    },
  ];
};

/**
 * takes in a string representation of a GenBank file and outputs our
 * part representation of it. an example of a Genbank file can be found
 * at ./parsers/Gebank, though there is significant variability to the
 * format
 *
 * another official example can be found at:
 * https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html
 */
export const genBankParser = (fileInput: string, fileName: string) =>
  fileInput
    .split(/\/\/\s/g)
    .filter((f) => f.length > 5)
    .map((file) => {
      // the first row contains the name of the part and its creation date
      // LOCUS       SCU49845     5028 bp    DNA             PLN       21-JUN-1999
      const HEADER_ROW = file.substring(
        file.indexOf("LOCUS"),
        file.search(/\\n|\n/)
      );
      const [, name] = HEADER_ROW.split(/\s{2,}/g).filter((h) => h);

      // trying to avoid giving a stupid name like Exported which Snapgene has by default
      // also, if there is not name in header, the seq length will be used as name, which should
      // be corrected (Number.parseInt to check for this case) https://stackoverflow.com/a/175787/7541747
      let parsedName = name;
      if (
        (parsedName === "Exported" && file.includes("SnapGene")) || // stupid Snapgene name
        Number.parseInt(parsedName, 10) // it thinks seq-length is the name
      ) {
        // first try and get the name from ACCESSION
        let accessionName = false;
        if (file.includes("ACCESSION")) {
          // this will be undefined is there is no
          const accession = file
            .substring(
              file.indexOf("ACCESSION"),
              file.indexOf("\n", file.indexOf("ACCESSION"))
            )
            .replace(".", "")
            .split(/\s{2,}/)
            .filter((a) => a !== "ACCESSION")
            .pop();
          if (accession) {
            parsedName = accession;
            accessionName = true;
          }
        }

        // otherwise, revert to trying to get the part name from the file name
        if (!accessionName && fileName) {
          parsedName = fileName
            .substring(
              0,
              Math.max(fileName.search(/\n|\||\./), fileName.lastIndexOf("."))
            )
            .replace(/\/\s/g, "");
        } else if (!accessionName) {
          parsedName = "Unnamed"; // give up
        }
      }

      // the part sequence is contained in and after the line that begins with ORIGIN
      // do this before annotations so we can calc seqlength
      //
      // ORIGIN
      //    1 gatcctccat atacaacggt atctccacct caggtttaga tctcaacaac ggaaccattg
      //    61 ccgacatgag acagttaggt atcgtcgaga gttacaagct aaaacgagca gtagtcagct
      const SEQ_ROWS = file.substring(
        file.lastIndexOf("ORIGIN") + "ORIGIN".length,
        file.length
      );
      let seq = SEQ_ROWS.replace(/[^gatc]/gi, "");
      ({ seq } = complement(seq)); // seq and compSeq

      return {
        name: parsedName.trim() || fileName,
        seq: seq,
        type: guessType(seq),
      };
    });
