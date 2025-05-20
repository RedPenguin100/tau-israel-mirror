"use client";
import { useCallback, useRef, useState } from "react";
import styles from "./uploader.module.css";

function FileUploader({ setFile, onReset, clearFileAfterUpload = false }) {
  const [isHovered, setIsHovered] = useState(false);
  const [fileName, setFileName] = useState(undefined);
  const dragArea = useRef();
  const fileInput = useRef();

  const onDragOver = useCallback((event) => {
    setIsHovered(true);
    dragArea?.current.classList.add(styles.active);
  });

  const onDragLeave = useCallback((event) => {
    setIsHovered(false);
    dragArea?.current.classList.remove(styles.active);
  });

  const onDropFile = useCallback((event) => {
    setIsHovered(false);
    dragArea?.current.classList.remove(styles.active);

    if (event?.dataTransfer?.files[0]) {
      const fileObj = event.dataTransfer.files[0];
      setFileName(fileObj?.name);
    }
  });

  const onBrowse = useCallback((event) => {
    event.preventDefault();
    fileInput.current.click();
  });

  const onChange = useCallback((event) => {
    const files = event.target?.files;
    const selectedFile = files[0];
    if (selectedFile) {
      setFileName(selectedFile?.name);
      setFile(selectedFile);
      if (clearFileAfterUpload) {
        fileInput.current.value = "";
        setFileName("");
      }
    }
  });

  const onResetFileInput = useCallback(() => {
    const resetValue = "";
    fileInput.current.value = resetValue;
    setFileName(resetValue);
    setFile(undefined);
    onReset();
  });

  return (
    <div
      ref={dragArea}
      className={styles.drag_area}
      onDrop={onDropFile}
      onDragOver={onDragOver}
      onDragLeave={onDragLeave}
    >
      <p>{!isHovered ? "Drag & Drop GENBANK \\ FASTA file" : "Release file"}</p>
      <span>or</span>
      <button onClick={onBrowse}>Browse Files</button>
      <p className={styles.support}>
        Supports: FASTA, FNA, FFN, FAA, FRN, FA, GENBANK, GBK, GB
      </p>
      <input
        ref={fileInput}
        type="file"
        name="file"
        id="file"
        accept="text/x-fasta, .fasta, .fna, .ffn, .faa, .frn, .fa, .genbank, .gbk, .gb"
        className={styles.file_input}
        onChange={onChange}
      />
      {fileName && (
        <div className={styles.selected_file_container}>
          <div className={styles.selected_file}>
            <i className="fa-solid fa-file-lines" />
            <p>{fileName}</p>
            <span
              className={styles.clear_file_input}
              onClick={onResetFileInput}
            >
              <i className="fa-solid fa-xmark" />
            </span>
          </div>
        </div>
      )}
    </div>
  );
}

export default FileUploader;
