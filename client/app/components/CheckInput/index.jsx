"use client";
import styles from "./check-input.module.scss";

function CheckInput({
  // type,
  name,
  id,
  label,
  checked,
  onChange = () => {},
}) {
  return (
    <li className={styles.row}>
      <input
        type="checkbox"
        name={name}
        id={id}
        className={styles.checkbox}
        checked={checked}
        onChange={onChange}
      ></input>
      <label htmlFor={id}>{label && label}</label>
    </li>
  );
}

export default CheckInput;
