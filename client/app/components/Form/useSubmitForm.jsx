import { useCallback } from "react";
import { setIsLoading } from "../../store";
import { BACKEND_URL } from "../../constants"; // ✅ using dynamic backend URL
import { getFastaFormat } from "../../utils";

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
  setAsoSequences
) => {
  return useCallback(
    async (event) => {
      event.preventDefault();

      setIsLoading(true);
      setDownloadFile(undefined);

      if (!isFormValid()) {
        setIsLoading(false);
        return;
      }

      try {
        // ✅ EDIT: Use dynamic backend URL instead of localhost
        const response = await fetch(`${BACKEND_URL}/gen_aso`, {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
          },
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
          setIsLoading(false);
          return;
        }

        const data = await response.json();

        if (!data || !data.asoSequence) {
          setErrors(["No ASO sequence returned from server"]);
          setIsLoading(false);
          return;
        }

        setAsoSequences(data.asoSequence || []);
        const asoFasta = getFastaFormat(data.asoSequence);
        const file = new Blob([asoFasta], { type: "text/plain" });
        setDownloadFile(file);

        setIsLoading(false);
      } catch (error) {
        setErrors([`Request failed: ${error.message}`]);
        setIsLoading(false);
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
    ]
  );
};
