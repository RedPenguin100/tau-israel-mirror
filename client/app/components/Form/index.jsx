"use client";
import { useCallback, useEffect, useState } from "react";
import FastaUploader from "../FastaUploader";
import styles from "./form.module.scss";
import Tabs from "../Tabs";
import CheckInput from "../CheckInput";
import RadioInput from "../RadioInput";

import ORGANISM from "../../data/organism";
import PLASMIDS from "../../data/plasmids";
import { isLoading } from "../../store";
import { useStore } from "@nanostores/react";
import { useSubmitForm } from "./useSubmitForm";

const INPUT_TYPES_ORGANISM = ["Select Organism", "Upload Organism Fasta"];
const INPUT_TYPES_GENE = ["Upload Gene Fasta", "Enter Gene Sequence"];


function isValidSequenceFile(filename) {
  const extension = filename.split('.').pop()?.toLowerCase();
  return [
    'fasta', 'fna', 'ffn', 'faa', 'frn',
    'fa', 'genbank', 'gbk', 'gb'
  ].includes(extension ?? '');
}


function Form({setAsoSequences}) {
  const [step, setStep] = useState(1);

  const [organismInputType, setOrganismInputType] = useState(INPUT_TYPES_ORGANISM[0]);
  const [geneInputType, setGeneInputType] = useState(INPUT_TYPES_GENE[0]);

  const [selectedOrganism, setSelectedOrganism] = useState(ORGANISM[0].id);
  const [organismFile, setOrganismFile] = useState(ORGANISM[0].file_name || "");
  const [geneFile, setGeneFile] = useState("none.fa"); // Default to no file
  const [geneSequence, setGeneSequence] = useState("ATGCTGACGTAGCTAGCTAGC");
  // Default gene sequence, can be replaced by user input

  const [numericParams, setNumericParams] = useState({ ASO_volume: 50, period_of_treatment: 7 });
  const [viewASO, setViewASO] = useState(true);
  const [errors, setErrors] = useState([]);
  const [isValid, setIsValid] = useState(true);
  const [downloadFile, setDownloadFile] = useState();
  
  const $isLoading = useStore(isLoading);

  const isFormValid = useCallback(() => {
    if (!organismFile.length) return setErrors(["Please provide organism sequence"]), false;
    if (!geneSequence.length) return setErrors(["Please provide gene sequence"]), false;

    const validNucleotides = ["A", "C", "T", "G"];
    if(!isValidSequenceFile(organismFile) ){
      setErrors(["Invalid organism sequence file format. Please upload a valid file."]);
      return false;
    }
      for (const char of new Set(geneSequence.toUpperCase())) {
        if (!validNucleotides.includes(char)) {
          setErrors(["Sequences must only contain nucleotides: A, C, G, T"]);
          return false;
        }
      }
    
    for (const [key, val] of Object.entries(numericParams)) {
      if (!(typeof val === "number" && val > 0)) {
        setErrors([`${key} must be a positive number`]);
        return false;
      }
    }

    setErrors([]);
    return true;
  }, [organismFile, geneSequence, numericParams]);

  useEffect(() => {
    setIsValid(isFormValid());
  }, [organismFile, geneSequence, numericParams, isFormValid]);

  const handleNumericParamChange = (paramName, value) => {
    setNumericParams((prev) => ({ ...prev, [paramName]: Number(value) }));
  };

  const submitForm = useSubmitForm(
    geneInputType,
    organismFile,
    geneFile,
    geneSequence,
    numericParams,
    viewASO,
    setDownloadFile,
    isFormValid,
    setErrors,
    setAsoSequences
  );

  return (
    <form onSubmit={submitForm} className={styles.form}>
      {step === 1 && (
        <>

          <h2>Organism Input</h2> 
          <Tabs
            options={INPUT_TYPES_ORGANISM}
            selectedOption={organismInputType}
            setSelectedOption={setOrganismInputType}
          />
          <section className={styles.sequence_input_1}>
            <div
              className={`${styles.container} ${organismInputType !== INPUT_TYPES_ORGANISM[0] && styles.disabled}`}
              onClick={() => setOrganismInputType(INPUT_TYPES_ORGANISM[0])}
            >
              <RadioInput
                name="selected-organism"
                options={ORGANISM}
                onChange={(selected) => {
                  setSelectedOrganism(selected);
                  const organism_data = ORGANISM.find((p) => p.id === selected);
                  setOrganismFile(organism_data?.file_name || "");
                }}
                defaultOption={selectedOrganism}
              />
            </div>
            <div
              className={`${styles.container} ${organismInputType !== INPUT_TYPES_ORGANISM[1] && styles.disabled}`}
              onClick={() => setOrganismInputType(INPUT_TYPES_ORGANISM[1])}
            >
              <FastaUploader setName={setOrganismFile} setSequence={()=>{}} clearFileAfterUpload />
            </div>
          </section>
        </>
      )}

      {step === 2 && (
        <>
          <h2>Gene Input</h2> 
          <Tabs
            options={INPUT_TYPES_GENE}
            selectedOption={geneInputType}
            setSelectedOption={setGeneInputType}
          />
          <section className={styles.sequence_input_2}>
            <div
              className={`${styles.container} ${geneInputType !== INPUT_TYPES_GENE[0] && styles.disabled}`}
              onClick={() => setGeneInputType(INPUT_TYPES_GENE[0])}
            >
              <FastaUploader setName={setGeneFile} setSequence={setGeneSequence} clearFileAfterUpload />
            </div>
            <div
              className={`${styles.container} ${geneInputType !== INPUT_TYPES_GENE[1] && styles.disabled}`}
              onClick={() => setGeneInputType(INPUT_TYPES_GENE[1])}
            >
              <textarea
                className={styles.text}
                placeholder="ACTG..."
                value={geneSequence}
                onChange={(e) => setGeneSequence(e.target.value)}
              />
            </div>
          </section>
        </>
      )}

      {step === 3 && (
        <>
          <h2>Numeric Parameters</h2> 
          <section className={styles.numeric_params}>
            {Object.entries(numericParams).map(([param, val]) => (
              <label key={param}>
                {param.replace(/_/g, " ")}:
                <input
                  type="number"
                  min={0}
                  value={val}
                  onChange={(e) => handleNumericParamChange(param, e.target.value)}
                />
              </label>
            ))}
          </section>

          <CheckInput
            type="checkbox"
            id="view-aso"
            name="view-aso"
            label="View ASO"
            checked={viewASO}
            onChange={(e) => setViewASO(e.target.checked)}
          />
        </>
      )}

      {errors.length > 0 && (
        <ul className={styles.errors}>
          {errors.map((error, idx) => (
            <li key={idx}>{error}</li>
          ))}
        </ul>
      )}

      <div className={styles.buttons}>
        {step > 1 && (
          <button
            type="button"
            className={styles.btn}
            onClick={() => setStep((s) => s - 1)}
          >
            Back
          </button>
        )}
        {step < 3 && (
          <button
            type="button"
            className={styles.btn}
            onClick={() => setStep((s) => s + 1)}
          >
            Next
          </button>
        )}
        {step === 3 && (
          <button
            className={`${styles.btn} ${styles.run_btn}`}
            disabled={!isValid || $isLoading}
            type="submit"
          >
            <span>Get Optimized ASO</span>
            <i className="fa-solid fa-arrow-right fa-lg" />
          </button>
        )}
        {downloadFile && (
          <a
            href={URL.createObjectURL(downloadFile)}
            download="aso_sequence.fa"
          >
            <button type="button" className={styles.btn}>
              Download ASO Sequence
              <i className="fa-solid fa-file-arrow-down fa-lg" />
            </button>
          </a>
        )}
      </div>
    </form>
  );
}

export default Form;
