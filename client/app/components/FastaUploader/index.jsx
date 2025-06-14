"use client";
import { useEffect, useState } from "react";
import FileUploader from "../FileUploader";

const CHUNK_SIZE = 1024 * 1024; // 1MB chunks

function FastaUploader({
  setName,
  setSequence,
  clearFileAfterUpload = false,
  onUpload = () => {},
  onReset = () => {},
}) {
  const [file, setFile] = useState();

  useEffect(() => {
    if (!file) return;

    const name = file.name.toLowerCase();
    if (name.endsWith(".fasta") || name.endsWith(".fa")|| name.endsWith(".fna")) {
      fastaMode(file);
    } else if (name.endsWith(".gb") || name.endsWith(".gbk")) {
      gbMode(file);
    } else {
      alert("Unsupported file type. Upload .fasta / .fa / .gb / .gbk");
    }
  }, [file]);

  // Helper to upload chunk with form data
  const uploadChunk = async (chunk, filename, isLastChunk, url) => {
    const formData = new FormData();
    formData.append("filename", filename);
    formData.append("is_last_chunk", isLastChunk ? "true" : "false");
    formData.append("chunk", chunk);

    const res = await fetch(url, {
      method: "POST",
      body: formData,
    });

    if (!res.ok) {
      const errorData = await res.json();
      throw new Error(errorData.error || "Upload failed");
    }
  };

  const fastaMode = async (file) => {
    // Use original file name + timestamp as unique filename on server
    const filename = `${file.name}_${Date.now()}.fasta`;
    try {
      let offset = 0;
      while (offset < file.size) {
        const chunk = file.slice(offset, offset + CHUNK_SIZE);
        offset += CHUNK_SIZE;
        const isLastChunk = offset >= file.size;

        await uploadChunk(chunk, filename, isLastChunk, "http://localhost:8000/upload_gene_fa");
      }
      setName(filename); 
      onUpload();
      if (clearFileAfterUpload) setFile(null);
    } catch (e) {
      alert("Upload error: " + e.message);
    }
  };

  const gbMode = async (file) => {
    const filename = `${file.name}_${Date.now()}.gb`;
    try {
      let offset = 0;
      while (offset < file.size) {
        const chunk = file.slice(offset, offset + CHUNK_SIZE);
        offset += CHUNK_SIZE;
        const isLastChunk = offset >= file.size;

        await uploadChunk(chunk, filename, isLastChunk, "http://localhost:8000/upload_gene_gb");
      }
      setName(filename); 
      onUpload();
      if (clearFileAfterUpload) setFile(null);
    } catch (e) {
      alert("Upload error: " + e.message);
    }
  };

  return (
    <FileUploader
      setFile={setFile}
      onReset={onReset}
      clearFileAfterUpload={clearFileAfterUpload}
    />
  );
}

export default FastaUploader;
