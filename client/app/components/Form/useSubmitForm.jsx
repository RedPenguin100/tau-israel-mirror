import { useCallback } from "react";
import { setIsLoading } from "../../store";
import { BACKEND_URL, EMAIL_BACKEND_URL } from "../../constants";
import { getFastaFormat } from "../../utils";

// helper: convert Blob to base64 string
async function blobToBase64(blob) {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onloadend = () => resolve(reader.result.split(",")[1]); // strip data: prefix
    reader.onerror = reject;
    reader.readAsDataURL(blob);
  });
}

export const useSubmitForm = (
  geneInputType,
  organismFile,
  geneFile,
  geneSequence,
  numericParams,
  viewASO,
  setDownloadFile,
  isFormValid,
  setErrors,
  setAsoSequences,
  userInfo,             // new for emails
  setShowThankYou,
  setIsProcessing
) => {
  return useCallback(
    async (event) => {
      event.preventDefault();

//      setIsLoading(true);
      setDownloadFile(undefined);

      if (!isFormValid()) {
//        setIsLoading(false);
        return;
      }

      try {
        // 1. Send "processing started" email WITH asoData (no result yet)
        setIsProcessing(true);
        await fetch(`${EMAIL_BACKEND_URL}/send-email`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            email: userInfo.email,
            name: userInfo.name,
            type: "processing_started",
            asoData: {
              geneInputType,
              organismFile,
              geneFile,
              geneSequence,
              numericParams,
              viewASO,
            },
          }),
        });
        setIsProcessing(false);
        setShowThankYou(true);
        // 2. Run gen_aso
        const response = await fetch(`${BACKEND_URL}/gen_aso`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            geneInputType,
            organismFile,
            geneFile,
            geneSequence,
            numericParams,
            viewASO,
          }),
        });

        if (!response.ok) {
          const errorText = await response.text();
          setErrors([`Server error: ${errorText}`]);
//          setIsLoading(false);
          setIsProcessing(false);
          return;
        }

        const data = await response.json();

        if (!data || !data.asoSequence) {
          setErrors(["No ASO sequence returned from server"]);
//          setIsLoading(false);
          setIsProcessing(false);
          return;
        }

        setAsoSequences(data.asoSequence || []);
        const asoFasta = getFastaFormat(data.asoSequence);
        const file = new Blob([asoFasta], { type: "text/plain" });
        setDownloadFile(file);

        // convert Blob -> base64
        const fileBase64 = await blobToBase64(file);

        // 3. Send "processing completed" email WITH file attached
        await fetch(`${EMAIL_BACKEND_URL}/send-email`, {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            email: userInfo.email,
            name: userInfo.name,
            type: "processing_completed",
            asoData: {
              geneInputType,
              organismFile,
              geneFile,
              geneSequence,
              numericParams,
              viewASO,
            },
            file: fileBase64, // base64 string
          }),
        });

//        setIsLoading(false);
      } catch (error) {
        console.error("Error:", error);
        setErrors([`Request failed: ${error.message}`]);
//       setIsLoading(false);
        setIsProcessing(false);
      }
    },
    [
      geneInputType,
      organismFile,
      geneFile,
      geneSequence,
      numericParams,
      viewASO,
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
