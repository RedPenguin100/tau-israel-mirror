import React from "react";
import styles from "./loading.module.css";

function LoadingPlasmid() {
  return (
    <div className={styles.modal}>
      <img
        className={styles.spinner}
        src="https://static.igem.wiki/teams/4648/wiki/software/plasmid.svg"
        width={100}
        height={100}
      />
    </div>
  );
}

export default LoadingPlasmid;
