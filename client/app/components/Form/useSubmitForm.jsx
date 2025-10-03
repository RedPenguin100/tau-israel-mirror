import { useCallback } from "react";
import { BACKEND_URL } from "../../constants";
import { getFastaFormat } from "../../utils";


export const useSubmitForm = (
  selectedGene,
  organismFile,
  geneSequence,
  topK,
  includeFeatureBreakdown,
  setDownloadFile,
  isFormValid,
  setErrors,
  setAsoSequences,
  userInfo,
  setShowThankYou,
  setIsProcessing
) => {
  return useCallback(
    async (event) => {
      event.preventDefault();
      setDownloadFile(undefined);

      if (!isFormValid()) return;

      setIsProcessing(true);

      try {
        // Single request to backend to handle everything
        const response = await fetch(`${BACKEND_URL}/run_aso`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            geneName: selectedGene,
            organismFile,
            geneSequence,
            top_k: topK,
            includeFeatureBreakdown,
            userEmail: userInfo.email,
            userName: userInfo.name,
          }),
        });

        if (!response.ok) {
          const errorText = await response.text();
          setErrors([`Server error: ${errorText}`]);
          setIsProcessing(false);
          return;
        }

        const data = await response.json();

        // Frontend gets minimal data (ASO sequences for preview / download)
        if (data?.asoSequence) {
          setAsoSequences(data.asoSequence);
          const asoFasta = getFastaFormat(data.asoSequence);
          const file = new Blob([asoFasta], { type: "text/plain" });
          setDownloadFile(file);
        }

        // Show thank you immediately
        setShowThankYou(true);
        setIsProcessing(false);
      } catch (error) {
        console.error("Error:", error);
        setErrors([`Request failed: ${error.message}`]);
        setIsProcessing(false);
      }
    },
    [
      selectedGene,
      organismFile,
      geneSequence,
      topK,
      includeFeatureBreakdown,
      setDownloadFile,
      isFormValid,
      setErrors,
      setAsoSequences,
      userInfo,
      setShowThankYou,
      setIsProcessing,
    ]
  );
};
