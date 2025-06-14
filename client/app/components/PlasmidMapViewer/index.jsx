"use client";
import { useEffect, useRef } from "react";
import { Circular, Linear, SeqViz } from "seqviz";
import styles from "./map.module.css";

//no need to edit

const DESKTOP_HEIGHT = "75vh";
const MOBILE_HEIGHT = "95vh";

function PlasmidMapViewer({ name, sequence, annotations, enzymes = [] }) {
  const circularRef = useRef();
  const linearRef = useRef();

  // scroll to bottom on mount
  useEffect(() => {
    circularRef.current?.scrollIntoView({ behavior: "smooth" });
  }, [circularRef]);

  return (
    <div className="plasmid-container">
    <SeqViz
      name={name}
      seq={sequence}
      annotations={annotations}
      style={{
        height: window.innerWidth < 600 ? MOBILE_HEIGHT : DESKTOP_HEIGHT,
        width: "100%",
      }}
      refs={{ circular: circularRef, linear: linearRef }}
    >
      {({ circularProps, linearProps, ...props }) => (
        <div className={styles.container}>
          <div ref={linearRef} className={styles.linear}>
            <Linear {...linearProps} {...props} />
          </div>
          <div ref={circularRef} className={styles.circle}>
            <Circular {...circularProps} {...props} />
          </div>
        </div>
      )}
    </SeqViz>
     </div>
  );
}

export default PlasmidMapViewer;
