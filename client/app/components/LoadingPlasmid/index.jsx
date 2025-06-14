import React, { useEffect, useState } from "react";
import styles from "./loading.module.css";
import "./DnaLoader.css";


//TO DO: find new logo for loading

function LoadingPlasmid() {
  const [progress, setProgress] = useState(0);
  const [startTime] = useState(Date.now());

  useEffect(() => {
    // Simulate backend load with random duration
    const duration = Math.random() * 3000 + 2000; // 2-5 seconds
    const interval = setInterval(() => {
      const elapsed = Date.now() - startTime;
      const percent = Math.min((elapsed / duration) * 100, 100);
      setProgress(percent);
      if (percent >= 100) clearInterval(interval);
    }, 50);
    return () => clearInterval(interval);
  }, [startTime]);

  return (
    <div className="dna-loader-wrapper">
      <div className="dna-helix">
        {[...Array(12)].map((_, i) => (
          <div key={i} className={`strand ${i % 2 === 0 ? "even" : "odd"}`} />
        ))}
      </div>
      <div className="progress-bar">
        <div className="progress-fill" style={{ width: `${progress}%` }} />
      </div>
      <div className="time-text">{progress.toFixed(0)}%</div>
    </div>
  );
}

export default LoadingPlasmid;
