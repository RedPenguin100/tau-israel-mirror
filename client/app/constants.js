// constants.js
const IS_LOCAL = process.env.NODE_ENV === "development";

export const BACKEND_URL = process.env.BACKEND_URL || (IS_LOCAL ? "http://localhost:8000" : "");
