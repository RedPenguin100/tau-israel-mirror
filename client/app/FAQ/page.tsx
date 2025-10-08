"use client";
import { useEffect, useState } from "react";
import Header from "../components/Header";

const faqData = [
  {
    question: "What is ASO Designer?",
    answer:
      "ASO Designer is a user-friendly tool that helps you design antisense oligonucleotides (ASOs) to silence specific genes. It uses advanced algorithms to suggest sequences that are efficient, specific, and practical for research use.",
  },
  {
    question: "Can I use ASO Designer for a different organism?",
    answer:
      "Right now, ASO Designer supports human, mouse, yeast and E. coli. These are the best-validated organisms in the tool. If you’re working with another species, feel free to send us an inquiry email at igem@tauex.tau.ac.il if you’d like us to add your organism of interest in future updates.",
  },
  {
    question: "How long does it usually take to get results?",
    answer:
      "For most genes, results appear within a few minutes up to half an hour. The exact time depends on the length of the gene, organism, the complexity of the structure, and the server availability.",
  },
  {
    question: "What factors determine the best ASOs that are being generated?",
    answer:
      "We rank ASOs based on a combination of: Binding energy (how strongly it pairs with the target), RNA folding (accessibility of the target site), GC content (stability balance), Off-target predictions and more features. This ensures that top suggestions are both effective and practical.",
  },
  {
    question: "Does ASO Designer guarantee experimental success?",
    answer:
      "No computational tool can fully guarantee success. ASO Designer provides high-quality predictions to save you time and resources, but all candidates should still be validated experimentally.",
  },
  {
    question: "What if my gene is not found in the database?",
    answer:
      "You can upload your own FASTA sequence, and ASO Designer will design ASOs directly from that input. You don’t need to rely only on our built-in databases.",
  },
  {
    question: "Can I design ASOs for non-coding RNAs or lncRNAs?",
    answer:
      "Yes. ASO Designer works with any transcript sequence, including coding and non-coding RNAs.",
  },
  {
    question: "Do I need advanced bioinformatics skills to use ASO Designer?",
    answer:
      "Not at all. The interface is designed to be intuitive for researchers, clinicians, and students. You don’t need coding skills, just the gene or sequence of interest.",
  },
  {
    question: "Is ASO Designer free to use?",
    answer:
      "Yes, ASO Designer is freely available for academic and research use. Please cite us in your publications if you use it in your work.",
  },
  {
    question: "Who do I contact if I have questions or need help?",
    answer:
      "We value every user. If you encounter issues, please reach out via igem@tauex.tau.ac.il. Your feedback helps us improve and ensures your research runs smoothly.",
  },
];

export default function FAQPage() {
  const [openIndex, setOpenIndex] = useState<number | null>(null);
  const [reloaded, setReloaded] = useState(false);

  useEffect(() => {
    // Only reload once to avoid infinite loop
    if (!reloaded) {
      setReloaded(true);
      window.location.reload();
    }
  }, [reloaded]);

  return (
    <div>
      <Header title="FAQ" subtitle="Frequently Asked Questions about ASO Designer" />
      <main style={{ padding: "1rem", maxWidth: "800px", margin: "0 auto" }}>
        {faqData.map((item, index) => (
          <div key={index} style={{ marginBottom: "1rem" }}>
            <h3
              onClick={() => setOpenIndex(openIndex === index ? null : index)}
              style={{
                cursor: "pointer",
                backgroundColor: "#6d81f0ff", 
                padding: "0.5rem 1rem",
                borderRadius: "5px",
              }}
            >
              {item.question}
            </h3>
            {openIndex === index && (
              <p
                style={{
                  padding: "0.5rem 1rem",
                  margin: 0,
                  backgroundColor: "#e6f0fb", // lighter blue for answer
                  borderLeft: "4px solid #357ABD", // darker accent
                }}
              >
                {item.answer}
              </p>
            )}
          </div>
        ))}
      </main>
    </div>
  );
}
