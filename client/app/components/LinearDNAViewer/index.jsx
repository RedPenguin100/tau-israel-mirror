import React, { useState } from "react";
import "./LinearDNAViewer.css";
import PlasmidMapViewer from "../PlasmidMapViewer";
export function LinearDNAViewer({ sequences }) {
  const [selected, setSelected] = useState(sequences[0]);

  return (
<div >
  <select
    
    value={selected.name}
    onChange={(e) => {
      const selectedSeq = sequences.find((s) => s.name === e.target.value);
      setSelected(selectedSeq);
    }}
  >
    {sequences.map((s) => (
      <option key={s.name} value={s.name}>
        {s.name}
      </option>
    ))}
  </select>

  <div >
                <PlasmidMapViewer
              name={selected.name}
              sequence={selected.sequence}
              annotations={ []}
            />
  </div>
</div>

  );
}
