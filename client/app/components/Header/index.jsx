import Link from "next/link";
import "./header.css";

function Header({ title, subtitle }) {
  return (
    <header>
      <div className="hero">
        <section className="top">
          <img
            src={`${process.env.ASSET_PREFIX || ""}/our_logo.png`}
            className="logo"
          />
          <nav>
            <ul>
              <li>
                <Link
                  href="https://2023.igem.wiki/tau-israel/software/#Instructions"
                  target="_blank"
                >
                  Instructions
                </Link>
              </li>
              <li>
                <Link href="https://2023.igem.wiki/tau-israel" target="_blank">
                  Wiki
                </Link>
              </li>
            </ul>
          </nav>
        </section>
        <h1 className="title">{title}</h1>
        <p>{subtitle}</p>
      </div>
    </header>
  );
}

export default Header;
