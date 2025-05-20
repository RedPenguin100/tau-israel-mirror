"use client";
import styles from "./page.module.css";
import Header from "./components/Header";
import Form from "./components/Form";
import LoadingPlasmid from "./components/LoadingPlasmid";
import PlasmidMapViewer from "./components/PlasmidMapViewer";
import { isLoading, plasmidMap } from "./store";
import { useStore } from "@nanostores/react";

export default function Home() {
  const $isLoading = useStore(isLoading);
  const $plasmidMap = useStore(plasmidMap);

  return (
    <>
      <Header
        title="Get Your Desired Plasmid Copy Number"
        subtitle="Using our tool, you can achieve a diverse range of plasmid copy numbers tailored to meet your exact requirements."
      />
      <div className={styles.main_wrapper}>
        <main className={styles.main}>
          <Form />
          {$isLoading && <LoadingPlasmid />}
          {$plasmidMap?.sequence && !$isLoading && (
            <PlasmidMapViewer
              name={$plasmidMap.name}
              sequence={$plasmidMap.sequence}
              annotations={$plasmidMap.annotations}
            />
          )}
        </main>
      </div>
    </>
  );
}
