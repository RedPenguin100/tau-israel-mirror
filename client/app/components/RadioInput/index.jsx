import React, { useEffect, useState } from "react";
import CheckInput from "../CheckInput";
import styles from "./radio.module.css";

function RadioInput({ name, options, onChange, defaultOption }) {
  const [selectedOption, setSelectedOption] = useState(defaultOption);
  useEffect(() => {
    onChange(selectedOption);
  }, []);
  return (
    <ul className={styles.container}>
      {options.map((option) => (
        <CheckInput
          key={option.id}
          name={name}
          id={option.id}
          label={option.label}
          checked={option.id == selectedOption}
          onChange={(event) => {
            const selectedId = event.target.id;
            if (event.target.checked) {
              setSelectedOption(selectedId);
              onChange(selectedId);
            }
          }}
          className={styles.option}
        />
      ))}
    </ul>
  );
}

export default RadioInput;
