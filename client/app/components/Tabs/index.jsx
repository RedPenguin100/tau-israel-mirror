import styles from "./tabs.module.scss";

function Tabs({ options, selectedOption, setSelectedOption, stacked = false }) {
  return (
    <section
      className={`${styles.switch_container}${
        stacked ? " " + styles.stacked : ""
      }`}
    >
      {options.map((option, index) => (
        <button
          key={index}
          type="button"
          onClick={(e) => setSelectedOption(e.target.value)}
          className={
            styles.switch_button +
            (option === selectedOption ? ` ${styles.active}` : "")
          }
          value={option}
        >
          {option}
        </button>
      ))}
      <input
        type="text"
        name="input-type"
        defaultValue={selectedOption}
        className={styles.input_type}
      />
    </section>
  );
}

export default Tabs;
