"use client";
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import FastaUploader from "../FastaUploader";
import styles from "./form.module.scss";
import CheckInput from "../CheckInput";
import RadioInput from "../RadioInput";
import { BACKEND_URL } from "../../constants";

import ORGANISM from "../../data/organism";
import { isLoading } from "../../store";
import { useStore } from "@nanostores/react";
import { useSubmitForm } from "./useSubmitForm";

function isValidSequenceFile(filename) {
  const extension = filename.split('.').pop()?.toLowerCase();
  return [
    'fasta', 'fna', 'ffn', 'faa', 'frn',
    'fa', 'genbank', 'gbk', 'gb'
  ].includes(extension ?? '');
}

function Form({ setAsoSequences }) {
  const [step, setStep] = useState(1);

  const [selectedOrganism, setSelectedOrganism] = useState(ORGANISM[0].id);
  const [organismFile, setOrganismFile] = useState(ORGANISM[0].file_name || "");

  const [customGeneSequence, setCustomGeneSequence] = useState("");
  const [customGeneFileName, setCustomGeneFileName] = useState("");
  const [geneNames, setGeneNames] = useState([]);
  const [geneSearchTerm, setGeneSearchTerm] = useState("");
  const [selectedGene, setSelectedGene] = useState("");
  const [geneNamesLoading, setGeneNamesLoading] = useState(false);
  const [geneNamesError, setGeneNamesError] = useState("");
  const [isFetchingSequence, setIsFetchingSequence] = useState(false);
  const [sequenceStatus, setSequenceStatus] = useState("");
  const sequenceStatusTimer = useRef(null);

  const [numericParams, setNumericParams] = useState({ ASO_volume: 50, period_of_treatment: 7 });
  const [viewASO, setViewASO] = useState(true);
  const [errors, setErrors] = useState([]);
  const [isValid, setIsValid] = useState(true);
  const [downloadFile, setDownloadFile] = useState();

  const [userInfo, setUserInfo] = useState({ name: "", email: "" });
  const [isProcessing, setIsProcessing] = useState(false);
  const [showThankYou, setShowThankYou] = useState(false);

  const $isLoading = useStore(isLoading);

  const isFormValid = useCallback((skipGeneValidation = false) => {
    if (!organismFile.length) return setErrors(["Please select an organism"]), false;
    if (!skipGeneValidation && !selectedGene) return setErrors(["Please select a gene from the list"]), false;

    if (!isValidSequenceFile(organismFile)) {
      setErrors(["Invalid organism sequence file format. Please select a valid file."]);
      return false;
    }

    for (const [key, val] of Object.entries(numericParams)) {
      if (!(typeof val === "number" && val > 0)) {
        setErrors([`${key} must be a positive number`]);
        return false;
      }
    }

    setErrors([]);
    return true;
  }, [organismFile, selectedGene, numericParams]);

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
      const skipGeneValidation = step < 2;
      setIsValid(isFormValid(skipGeneValidation));
    } else {
      setIsValid(isUserInfoValid());
    }
  }, [organismFile, numericParams, isFormValid, step, userInfo, isUserInfoValid, selectedGene]);

  useEffect(() => {
    let isMounted = true;

    const applyNames = (names) => {
      if (!isMounted) return;
      setGeneNames(names);
      setGeneNamesError(names.length ? "" : "No gene names available");
    };

    const fetchFromBackend = async () => {
      const response = await fetch(`${BACKEND_URL}/gene_names`);
      if (!response.ok) {
        throw new Error(`Backend responded with ${response.status}`);
      }
      const payload = await response.json();
      const names = Array.isArray(payload?.genes) ? payload.genes : [];
      applyNames(names);
    };

    const fetchFromStatic = async () => {
      const response = await fetch(`/gene_names.txt`);
      if (!response.ok) {
        throw new Error(`Static file responded with ${response.status}`);
      }
      const text = await response.text();
      const names = text
        .split(/\r?\n/)
        .map((line) => line.trim())
        .filter(Boolean);
      applyNames(names);
    };

    const fetchGeneNames = async () => {
      setGeneNamesLoading(true);
      try {
        await fetchFromBackend();
      } catch (backendError) {
        console.error("Failed to load gene names from backend", backendError);
        try {
          await fetchFromStatic();
        } catch (staticError) {
          console.error("Failed to load gene names from static fallback", staticError);
          if (isMounted) {
            setGeneNamesError("Failed to load gene list. Please refresh the page or try again later.");
          }
        }
      } finally {
        if (isMounted) {
          setGeneNamesLoading(false);
        }
      }
    };

    fetchGeneNames();

    return () => {
      isMounted = false;
    };
  }, []);

  const filteredGeneNames = useMemo(() => {
    if (!geneSearchTerm.trim()) return geneNames;
    const lowered = geneSearchTerm.trim().toLowerCase();
    return geneNames.filter((name) => name.toLowerCase().startsWith(lowered));
  }, [geneNames, geneSearchTerm]);

  useEffect(() => {
    return () => {
      if (sequenceStatusTimer.current) {
        clearTimeout(sequenceStatusTimer.current);
        sequenceStatusTimer.current = null;
      }
    };
  }, []);

  const handleNumericParamChange = (paramName, value) => {
    setNumericParams((prev) => ({ ...prev, [paramName]: Number(value) }));
  };

  const submitForm = useSubmitForm(
    selectedGene,
    organismFile,
    customGeneSequence,
    numericParams,
    viewASO,
    setDownloadFile,
    isFormValid,
    setErrors,
    setAsoSequences,
    userInfo,
    setShowThankYou,
    setIsProcessing
  );


  return (
    <form onSubmit={submitForm} className={styles.form}>
      {step === 1 && (
        <>
          <h2>Organism Input</h2> 
          <section className={styles.sequence_input_1}>
            <div className={`${styles.container} ${styles.wide}`}>
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
          </section>
        </>
      )}

      {step === 2 && (
        <>
          <h2>Gene Input</h2> 
          <section className={styles.sequence_input_2}>
            <div className={styles.container}>
              <div className={styles.geneListWrapper}>
                <input
                  type="text"
                  value={geneSearchTerm}
                  onChange={(event) => setGeneSearchTerm(event.target.value)}
                  placeholder="Search gene"
                  className={styles.geneSearch}
                />
                <div className={styles.geneList}>
                  {geneNamesLoading && <p>Loading gene names…</p>}
                  {!geneNamesLoading && geneNamesError && <p className={styles.geneListError}>{geneNamesError}</p>}
                  {!geneNamesLoading && !geneNamesError && filteredGeneNames.length === 0 && (
                    <p>No genes match the current search.</p>
                  )}
                  {!geneNamesLoading && !geneNamesError && filteredGeneNames.map((name) => (
                    <button
                      key={name}
                      type="button"
                      className={`${styles.geneListItem} ${selectedGene === name ? styles.geneListItemSelected : ""}`}
                      onClick={() => {
                        setSelectedGene(name);
                        setCustomGeneSequence("");
                        setCustomGeneFileName("");
                        setErrors([]);
                        if (sequenceStatusTimer.current) {
                          clearTimeout(sequenceStatusTimer.current);
                          sequenceStatusTimer.current = null;
                        }
                        setIsFetchingSequence(true);
                        setSequenceStatus("Loading sequence…");
                        sequenceStatusTimer.current = setTimeout(() => {
                          setIsFetchingSequence(false);
                          setSequenceStatus("Sequence loaded");
                          sequenceStatusTimer.current = null;
                        }, 800);
                      }}
                    >
                      {name}
                    </button>
                  ))}
                </div>
                {selectedGene && (
                  <div className={styles.selectedGene}>
                    Selected gene: <strong>{selectedGene}</strong>
                    {sequenceStatus && (
                      <span
                        className={`${styles.sequenceStatus} ${
                          isFetchingSequence
                            ? styles.sequenceStatusLoading
                            : customGeneSequence
                              ? styles.sequenceStatusCustom
                              : styles.sequenceStatusLoaded
                        }`}
                      >
                        {sequenceStatus}
                      </span>
                    )}
                  </div>
                )}
                {step === 2 && errors.length > 0 && (
                  <div className={styles.inlineErrors}>
                    {errors.map((error, index) => (
                      <span key={index} className={styles.inlineErrorItem}>{error}</span>
                    ))}
                  </div>
                )}
              </div>
            </div>
            <div className={styles.container}>
              <FastaUploader
                setName={setCustomGeneFileName}
                setSequence={(sequence) => {
                  if (sequenceStatusTimer.current) {
                    clearTimeout(sequenceStatusTimer.current);
                    sequenceStatusTimer.current = null;
                  }
                  setCustomGeneSequence(sequence);
                  setIsFetchingSequence(false);
                  setSequenceStatus(sequence ? "Custom sequence loaded" : "");
                  setErrors([]);
                }}
                clearFileAfterUpload
              />
              {(customGeneFileName || customGeneSequence) && (
                <div className={styles.selectedGene}>
                  Using uploaded sequence{customGeneFileName ? ` (${customGeneFileName})` : ""} – {customGeneSequence.length} nt
                  <button
                    type="button"
                    className={styles.clearOverride}
                    onClick={() => {
                      setCustomGeneSequence("");
                      setCustomGeneFileName("");
                      if (sequenceStatusTimer.current) {
                        clearTimeout(sequenceStatusTimer.current);
                        sequenceStatusTimer.current = null;
                      }
                      setIsFetchingSequence(false);
                      if (selectedGene) {
                        setSequenceStatus("Sequence loaded");
                      } else {
                        setSequenceStatus("");
                      }
                      setErrors([]);
                    }}
                  >
                    Use selected gene only
                  </button>
                </div>
              )}
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

      {errors.length > 0 && step !== 2 && (
        <ul className={styles.errors}>
          {errors.map((error, idx) => (
            <li key={idx}>{error}</li>
          ))}
        </ul>
      )}

      <div className={styles.buttons} style={{ display: 'flex', gap: '10px' }}>
        {step > 1 && !isProcessing && !showThankYou && (
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
            onClick={() => window.location.href = "/software-tools/tau-israel/"}
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
