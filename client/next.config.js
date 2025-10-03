/**
 * @type {import('next').NextConfig}
 */
const nextConfig = {
  reactStrictMode: true,
  images: {
    unoptimized: true,
  },
  output: "export",
  assetPrefix: process.env.ASSET_PREFIX || "",
  env: {
    NEXT_PUBLIC_BACKEND_URL: process.env.BACKEND_URL || "",
  },
};

module.exports = nextConfig;
