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
  
  // New state for user info
  const [userInfo, setUserInfo] = useState({ name: "", email: "" });
  const [isProcessing, setIsProcessing] = useState(false);
  const [showThankYou, setShowThankYou] = useState(false);
  
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

  const isUserInfoValid = useCallback(() => {
    if (!userInfo.name.trim()) return setErrors(["Please provide your name"]), false;
    if (!userInfo.email.trim()) return setErrors(["Please provide your email address"]), false;
    
    const emailRegex = /^[^\s@]+@[^\s@]+\.[^\s@]+$/;
    if (!emailRegex.test(userInfo.email)) {
      setErrors(["Please enter a valid email address"]);
      return false;
    }

    setErrors([]);
    return true;
  }, [userInfo]);

  useEffect(() => {
    if (step <= 3) {
      setIsValid(isFormValid());
    } else {
      setIsValid(isUserInfoValid());
    }
  }, [organismFile, geneSequence, numericParams, isFormValid, step, userInfo, isUserInfoValid]);

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

  const handleFormSubmit = async (e) => {
    e.preventDefault();
    
    if (step === 4 && isUserInfoValid()) {
      // Start processing and send emails
      setIsProcessing(true);
      
      try {
        // Send "processing started" email
        await fetch('/api/send-email', {
          method: 'POST',
          headers: {
            'Content-Type': 'application/json',
          },
          body: JSON.stringify({
            email: userInfo.email,
            name: userInfo.name,
            type: 'processing_started'
          }),
        });

        // Simulate processing time (in real app, this would be actual processing)
        setTimeout(async () => {
          // Send "processing completed" email
          await fetch('/api/send-email', {
            method: 'POST',
            headers: {
              'Content-Type': 'application/json',
            },
            body: JSON.stringify({
              email: userInfo.email,
              name: userInfo.name,
              type: 'processing_completed',
              asoData: {
                geneInputType,
                organismFile,
                geneFile,
                geneSequence,
                numericParams,
                viewASO
              }
            }),
          });
          
          setIsProcessing(false);
          setShowThankYou(true);
          // You could also call the actual ASO processing function here
        }, 5000);
        
      } catch (error) {
        console.error('Error sending emails:', error);
        setIsProcessing(false);
      }
    }
  };

  return (
    <form onSubmit={handleFormSubmit} className={styles.form}>
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

      {step === 4 && !showThankYou && (
        <>
          <h2>Contact Information</h2>
          <p>We'll send you an email when your ASO analysis is complete.</p>
          <section className={styles.user_info}>
            <div className={styles.formGroup}>
              <label htmlFor="name">Full Name *</label>
              <input
                type="text"
                id="name"
                name="name"
                value={userInfo.name}
                onChange={(e) => setUserInfo(prev => ({ ...prev, name: e.target.value }))}
                placeholder="Enter your full name"
                className={styles.userInput}
              />
            </div>

            <div className={styles.formGroup}>
              <label htmlFor="email">Email Address *</label>
              <input
                type="email"
                id="email"
                name="email"
                value={userInfo.email}
                onChange={(e) => setUserInfo(prev => ({ ...prev, email: e.target.value }))}
                placeholder="Enter your email address"
                className={styles.userInput}
              />
            </div>
          </section>
        </>
      )}

      {showThankYou && (
        <>
          <h2>Thank You!</h2>
          <div className={styles.thankYouMessage}>
            <div className={styles.successIcon}>
              <i className="fa-solid fa-check-circle fa-3x"></i>
            </div>
            <p>Thank you for your input, {userInfo.name}!</p>
            <p>Please make sure you received a confirmation email at <strong>{userInfo.email}</strong>.</p>
            <p>We'll send you another email with your ASO results once the analysis is complete.</p>
          </div>
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
        {step < 3 && (
          <button
            type="button"
            className={styles.btn}
            onClick={() => setStep((s) => s + 1)}
          >
            Next
          </button>
        )}
        {step > 1 && !isProcessing && (
          <button
            type="button"
            className={styles.btn}
            onClick={() => setStep((s) => s - 1)}
          >
            Back
          </button>
        )}
        {step === 3 && (
          <button
            type="button"
            className={`${styles.btn} ${styles.run_btn}`}
            disabled={!isValid || $isLoading}
            onClick={() => setStep(4)}
          >
            <span>Get Optimized ASO</span>
            <i className="fa-solid fa-arrow-right fa-lg" />
          </button>
        )}
        {step === 4 && !showThankYou && (
          <button
            className={`${styles.btn} ${styles.run_btn}`}
            disabled={!isValid || isProcessing}
            type="submit"
          >
            {isProcessing ? (
              <>
                <i className="fa-solid fa-spinner fa-spin" />
                <span>Processing...</span>
              </>
            ) : (
              <>
                <span>Start Processing</span>
                <i className="fa-solid fa-paper-plane fa-lg" />
              </>
            )}
          </button>
        )}

        {showThankYou && (
          <button
            type="button"
            className={styles.btn}
            onClick={() => window.location.href = "/"}
          >
            <i className="fa-solid fa-home fa-lg" />
            <span>Back to Home</span>
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
