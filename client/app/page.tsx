"use client";
import styles from "./page.module.css";
import Header from "./components/Header";
import Form from "./components/Form";

import LoadingPlasmid from "./components/LoadingPlasmid";
import {LinearDNAViewer} from "./components/LinearDNAViewer";
import { isLoading, plasmidMap } from "./store";
import { useStore } from "@nanostores/react";
import { useState } from "react";
export default function Home() {
  const $isLoading = useStore(isLoading);
  const $plasmidMap = useStore(plasmidMap);
  const[aso_sequences, setAsoSequences] = useState([]);
  

  return (
    <>
      <Header
        title="ASO Generator for Targeted Gene Knockdown"
        subtitle="Input a gene and organism to get a high-efficiency antisense oligo, instantly."
      />
      <div className={styles.main_wrapper}>
        <main className={styles.main}>
          <Form setAsoSequences={setAsoSequences}/>
          {$isLoading && <LoadingPlasmid />}
          { !$isLoading && aso_sequences.length > 0 &&
            (<LinearDNAViewer sequences = {aso_sequences} />) }


        </main>
      </div>
    </>
  );
}
