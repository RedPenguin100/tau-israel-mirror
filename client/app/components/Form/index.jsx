"use client";
import { useCallback, useEffect, useState } from "react";
import FastaUploader from "../FastaUploader";
import styles from "./form.module.scss";
import Tabs from "../Tabs";
import CheckInput from "../CheckInput";
import RadioInput from "../RadioInput";
import PLASMIDS from "../../data/plasmids";
import { isLoading } from "../../store";
import { useStore } from "@nanostores/react";
import { useSubmitForm } from "./useSubmitForm";
import { getSingleCopyNumber } from "@/app/utils";

const INPUT_TYPES = ["Selected Plasmid", "Sequence File", "Sequence as Text"];

function Form() {
  const [inputType, setInputType] = useState(INPUT_TYPES[0]);
  const [selectedPlasmid, setSelectedPlasmid] = useState(PLASMIDS[0].id);
  const [name, setName] = useState("unknown");
  const [sequence, setSequence] = useState("");
  const [copyNumberInputType, setCopyNumberInputType] = useState("single"); // 'single' or 'grid'
  const [copyNumber, setCopyNumber] = useState([50, 100, 300, 500, 800]);
  const [doVisualize, setDoVisualize] = useState(true);
  const [isValid, setIsValid] = useState(true);
  const [errors, setErrors] = useState([]);
  const [downloadFile, setDownloadFile] = useState();
  const $isLoading = useStore(isLoading);

  useEffect(() => {
    const tempIsValid = isFormValid();
    if (tempIsValid != isValid) {
      setIsValid(tempIsValid);
    }
  }, [inputType, name, sequence, copyNumber]);

  const isFormValid = useCallback(() => {
    if (sequence.length == 0) {
      setErrors(["please fill plasmid sequence"]);
      return false;
    }
    for (const i of new Set(sequence)) {
      if (!["A", "C", "T", "G"].includes(i.toUpperCase(i))) {
        setErrors([
          "sequence must be constructed of nucleotides only: A, C, G, T",
        ]);
        return false;
      }
    }
    if (!copyNumber.every((pcn) => pcn > 0 && pcn <= 900)) {
      setErrors(["desired copy number must be in the following range 1 - 900"]);
      return false;
    }
    setErrors([]);
    return true;
  }, [inputType, name, sequence, copyNumber, doVisualize]);

  const submitForm = useSubmitForm(
    inputType,
    name,
    sequence,
    copyNumber,
    copyNumberInputType,
    doVisualize,
    setDownloadFile,
    isFormValid,
    setErrors
  );

  return (
    <form onSubmit={submitForm} className={styles.form}>
      <span className={styles.sequence_input_tabs}>
        <Tabs
          options={INPUT_TYPES}
          selectedOption={inputType}
          setSelectedOption={setInputType}
        />
      </span>
      <section
        className={`${styles.sequence_input} ${
          inputType === INPUT_TYPES[0]
            ? styles.plasmid
            : inputType === INPUT_TYPES[1]
            ? styles.fasta
            : styles.sequence
        }`}
      >
        <div
          className={`${styles.container} ${
            inputType != INPUT_TYPES[0] && styles.disabled
          }`}
          onClick={() => setInputType(INPUT_TYPES[0])}
        >
          <RadioInput
            name="selected-plasmid"
            options={PLASMIDS}
            onChange={(selected) => {
              setSelectedPlasmid(selected);
              const plasmid_data = PLASMIDS.filter((p) => p.id == selected)[0];
              setName(plasmid_data.label);
              setSequence(plasmid_data.sequence);
            }}
            defaultOption={selectedPlasmid}
          />
        </div>
        <div
          className={`${styles.container} ${
            inputType != INPUT_TYPES[1] && styles.disabled
          }`}
          onClick={(event) => {
            if (inputType != INPUT_TYPES[1]) {
              event.preventDefault();
              setInputType(INPUT_TYPES[1]);
            }
          }}
        >
          <FastaUploader
            setName={setName}
            setSequence={setSequence}
            clearFileAfterUpload
            onUpload={() => setInputType(INPUT_TYPES[2])}
            onReset={() => {
              setName("");
              setSequence("");
            }}
          />
        </div>
        <div
          className={`${styles.container} ${
            inputType != INPUT_TYPES[2] && styles.disabled
          }`}
          onClick={() => setInputType(INPUT_TYPES[2])}
        >
          <textarea
            className={styles.text}
            placeholder="ACTG...|"
            value={sequence}
            onChange={(e) => setSequence(e.target.value)}
            id="sequence"
          />
        </div>
      </section>
      {errors.length > 0 && (
        <ul className={styles.errors}>
          {errors.map((error, index) => (
            <p key={index}>{error}</p>
          ))}
        </ul>
      )}

      <div className={styles.pcn_row}>
        <span className={styles.pcn_inputs_container}>
          <label htmlFor="pcn-0" className={styles.pcn_label}>
            Desired Copy Number:
          </label>
          <div className={styles.pcn_inputs}>
            {(copyNumberInputType == "single"
              ? getSingleCopyNumber(copyNumber)
              : copyNumber
            ).map((pcn, index) => (
              <input
                type="number"
                id={`pcn-${index}`}
                key={index}
                name="pcn"
                value={pcn}
                onChange={(e) => {
                  const updatedCopyNumber = [...copyNumber];
                  updatedCopyNumber[index] = e.target.value;
                  setCopyNumber(updatedCopyNumber);
                }}
                className={styles.pcn_input}
                max={900}
                min={0}
              />
            ))}
          </div>
        </span>
        <Tabs
          options={["single", "grid"]}
          selectedOption={copyNumberInputType}
          setSelectedOption={setCopyNumberInputType}
          stacked
        />
      </div>
      <div className={styles.row}>
        <CheckInput
          type="checkbox"
          id="visualize"
          name="visualize"
          label="visualize the plasmid map"
          checked={doVisualize}
          onChange={(event) => setDoVisualize(event.target.checked)}
        />
      </div>

      <div className={styles.buttons}>
        <button
          className={`${styles.btn} ${styles.run_btn}`}
          disabled={!isValid || $isLoading}
          type="submit"
        >
          <span>Get Modified Plasmid</span>
          <i className="fa-solid fa-arrow-right fa-lg" />
        </button>
        {downloadFile && (
          <a
            href={URL.createObjectURL(downloadFile)}
            download={`${name}_pcn.fa`}
          >
            <button type="button" className={styles.btn}>
              Download Modified Sequence
              <i className="fa-solid fa-file-arrow-down fa-lg" />
            </button>
          </a>
        )}
      </div>
    </form>
  );
}

export default Form;
