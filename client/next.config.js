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
};

module.exports = nextConfig;
