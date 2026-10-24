import { fileURLToPath } from "node:url";
export default {
  root: fileURLToPath(new URL(".", import.meta.url)),
  publicDir: fileURLToPath(new URL("./public", import.meta.url)),
  // no live reload: an edit must never navigate a page that is being rendered
  server: { port: 5241, strictPort: true, hmr: false, fs: { allow: [fileURLToPath(new URL("..", import.meta.url))] } },
};
