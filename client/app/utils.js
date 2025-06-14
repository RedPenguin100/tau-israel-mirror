export const getSingleCopyNumber = (copyNumber) => copyNumber.slice(0, 1);

export const replaceSubSequence = (sequence, subsequence, start, end) =>
  sequence.slice(0, start) + subsequence + sequence.slice(end);

export const handleModifications = (name, original_sequence, modifications) => {
  return modifications.map(({ subsequence, start, end, pcn }) => ({
    name: `${name}${pcn ? " copy number " + pcn : ""}`,
    sequence: replaceSubSequence(original_sequence, subsequence, start, end),
  }));
};

export const getFastaFormat = (sequences) => {
  return sequences
    .map(({ name, sequence }) => `>${name}\n${sequence}`)
    .join("\n");
};
