import { useCallback } from "react";
import { resetPlasmidMap, setIsLoading, setPlasmidMap } from "../../store";
import { BACKEND_URL } from "../../constants";
import {
  getFastaFormat,
  getSingleCopyNumber,
  handleModifications,
} from "../../utils";

export const useSubmitForm = (
  inputType,
  name,
  sequence,
  copyNumber,
  copyNumberInputType,
  doVisualize,
  setDownloadFile,
  isFormValid,
  setErrors
) => {
  return useCallback(
    async (event) => {
      event.preventDefault();

      setIsLoading(true);
      setDownloadFile(undefined);
      resetPlasmidMap();

      if (isFormValid()) {
        // https://developer.mozilla.org/en-US/docs/Web/API/Fetch_API/Using_Fetch
        try {
          const response = await fetch(BACKEND_URL, {
            method: "POST",
            headers: {
              "Content-Type": "application/json",
            },
            body: JSON.stringify({
              sequence,
              copy_number:
                copyNumberInputType == "single"
                  ? getSingleCopyNumber(copyNumber)
                  : copyNumber,
              visualize_plasmid_map: doVisualize,
            }),
          });

          const responseBody = await response.json();
          setErrors(responseBody.errors);
          if (
            !responseBody.errors.length &&
            responseBody.modifications.length
          ) {
            const modifiedSequences = handleModifications(
              name,
              sequence,
              responseBody.modifications
            );
            setDownloadFile(
              new File(getFastaFormat(modifiedSequences), "foo.txt", {
                type: "text/plain",
              })
            );
            if (doVisualize) {
              setPlasmidMap({
                name,
                sequence: modifiedSequences.at(-1).sequence,
                annotations: responseBody.annotations,
              });
            }
          }
        } catch (error) {
          console.log({ error });
          setErrors([error.message]);
        } finally {
          setIsLoading(false);
        }
      }
    },
    [inputType, name, sequence, copyNumber, copyNumberInputType, doVisualize]
  );
};
