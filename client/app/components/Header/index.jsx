import Link from "next/link";
import "./header.css";

function Header({ title, subtitle }) {
  return (
    <header>
        <section className="top">
          <img
            src={`${process.env.ASSET_PREFIX || ""}/our_logo.png`}
            className="logo"
          />
          <nav>
            <ul>
              <li>
                <Link
                  href="https://2025.igem.wiki/tau-israel/software"
                  target="_blank"
                  className="my-link"
                >
                  Software
                </Link>
              </li>
              <li>
                <Link href="https://2025.igem.wiki/tau-israel" 
                      target="_blank"
                      className="my-link"                  
                                        >
                  Wiki
                </Link>
              </li>
            </ul>
          </nav>
        </section>
        <h1 className="title">{title}</h1>
        <p>{subtitle}</p>
    </header>
  );
}

export default Header;
