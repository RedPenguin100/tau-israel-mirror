import Header from "../components/Header";
import styles from "./user_guide.module.css";

const steps = [
  {
    title: "1. Select an organism",
    details: [
      "Pick one of the optional organisms to preload the matching reference sequence.",
    ],
  },
  {
    title: "2. Provide the target gene",
    details: [
      "Upload a FASTA/GenBank file or paste the nucleotide sequence directly.",
      "The designer accepts only A, C, G, T characters, clean the sequence before pasting."
    ],
  },
  {
    title: "3. Tune numeric parameters",
    details: [
      "Set ASO volume and treatment period to guide downstream recommendations.",
      "Leave the defaults if you are evaluating the pipeline for the first time."
    ],
  },
  {
    title: "4. Enter contact details and run",
    details: [
      "Provide your name and email so the pipeline can notify you when processing is complete.",
      "Click \"Start Processing\" to queue the job. You will receive a confirmation email that your sequence is being processed."
    ],
  },
  {
    title: "5. Review your output",
    details: [
      "Once complete, you receive the result email with a downloadable FASTA file of the generated ASO candidate.",
      "The page also displays the sequences for quick inspection and visualization."
    ],
  },
];

const tips = [
  "If processing appears stuck, refresh the page and check your inbox for the confirmation email.",
];

export default function UserGuidePage() {
  return (
    <div>
      <Header title="User Guide" subtitle="How to use the ASO Designer tool" />
      <main className={styles.wrapper}>
        <section className={styles.intro}>
          <p>
            The TA(U)SO designer walks you through a short four-step form. Follow the steps below to
            submit a sequence and receive antisense oligonucleotide (ASO) candidates tailored to your
            organism of interest and target gene.
          </p>
        </section>

        <section className={styles.steps}>
          {steps.map(({ title, details }) => (
            <article key={title} className={styles.stepCard}>
              <h2>{title}</h2>
              <ul>
                {details.map((line) => (
                  <li key={line}>{line}</li>
                ))}
              </ul>
            </article>
          ))}
        </section>

        <section className={styles.tips}>
          <h2>Quick tip</h2>
          <ul>
            {tips.map((tip) => (
              <li key={tip}>{tip}</li>
            ))}
          </ul>
        </section>
      </main>
    </div>
  );
}
