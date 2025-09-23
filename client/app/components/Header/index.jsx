import Link from "next/link";
import Image from "next/image";
import ourLogo from "../public/our_logo.png"; // import from public folder
import "./header.css";

function Header({ title, subtitle }) {
  return (
    <header>
      {/* Wave layers, May addition */}
      <img src="https://static.igem.wiki/teams/5661/waves/wave-1.svg" alt="" className="waves"/>
      <img src="https://static.igem.wiki/teams/5661/waves/wave-2.svg" alt="" className="waves"/>
      <img src="https://static.igem.wiki/teams/5661/waves/wave-3.svg" alt="" className="waves"/>
      <img src="https://static.igem.wiki/teams/5661/waves/wave-4.svg" alt="" className="waves"/>
      <img src="https://static.igem.wiki/teams/5661/waves/wave-5.svg" alt="" id="shape" />

      <section className="top">
        <Image
          src={ourLogo}
          alt="Logo"
          className="logo"
          width={100} // adjust to your desired size
          height={100} // adjust to your desired size
        />
        <nav>
          <ul>
            <li>
              <Link href="/tau-israel/software/" className="my-link">ASO Designer</Link>
            </li>
            <li>
              <Link href="/tau-israel/software/FAQ" className="my-link">FAQ</Link>
            </li>
            <li>
              <Link href="/tau-israel/software/user_guide" className="my-link">User Guide</Link>
            </li>
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
              <Link
                href="https://2025.igem.wiki/tau-israel"
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
