import "./globals.css";
import type { Metadata } from "next";
import { Inter } from "next/font/google";

const inter = Inter({ subsets: ["latin"] });

export const metadata: Metadata = {
  title: "PCN Gradient Tool",
  description: "Plasmid Copy Number Control Tool",
};

export default function RootLayout({
  children,
}: {
  children: React.ReactNode;
}) {
  return (
    <html lang="en">
      <head>
        <link
          href={`${
            process.env.ASSET_PREFIX || ""
          }/fontawesome/css/fontawesome.css`}
          rel="stylesheet"
        />
        <link
          href={`${process.env.ASSET_PREFIX || ""}/fontawesome/css/brands.css`}
          rel="stylesheet"
        />
        <link
          href={`${process.env.ASSET_PREFIX || ""}/fontawesome/css/solid.css`}
          rel="stylesheet"
        />
      </head>
      <body className={inter.className}>{children}</body>
    </html>
  );
}
