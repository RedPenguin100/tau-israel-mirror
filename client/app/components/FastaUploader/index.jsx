"use client";
import { useEffect, useState } from "react";
import FileUploader from "../FileUploader";
import { fastaParser, genBankParser } from "../../geneFileParsers";

const fastaExtensions = ["fasta", "fna", "ffn", "faa", "frn", "fa"];
const geneBankExtensions = ["genbank", "gbk", "gb"];

function FastaUploader({
  setName,
  setSequence,
  clearFileAfterUpload = false,
  onUpload = () => {},
  onReset = () => {},
}) {
  const [file, setFile] = useState();

  useEffect(() => {
    if (file) {
      const reader = new FileReader();
      reader.onload = async (e) => {
        const fileContent = e.target.result;
        const fileExtension = file.name.split(".").at(-1).toLowerCase();
        let sequences;
        if (fastaExtensions.includes(fileExtension)) {
          sequences = fastaParser(fileContent, file.name);
        } else if (geneBankExtensions.includes(fileExtension)) {
          sequences = genBankParser(fileContent, file.name);
        } else {
          sequences = [
            {
              seq: `Something wrong has happened, you may have uploaded an unsupported file: .${fileExtension}`,
              name: "error",
            },
          ];
        }

        let { seq, name } = sequences[0];

        // clear sequence from wrong charachters
        if ([...fastaExtensions, ...geneBankExtensions].includes()) {
          seq = seq.replace(/[^ACGT]/gim, "").toLocaleUpperCase();
        }

        setName(name);
        setSequence(seq);
        onUpload();
      };

      reader.readAsText(file);
    }
  }, [file]);
  return (
    <FileUploader
      setFile={setFile}
      onReset={onReset}
      clearFileAfterUpload={clearFileAfterUpload}
    />
  );
}

export default FastaUploader;