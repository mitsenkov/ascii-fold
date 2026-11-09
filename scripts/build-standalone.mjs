// Builds dist/ascii-fold-standalone.html: the full app in one file (CSS + JS inlined),
// so it can be opened straight from disk (file://) or dropped on any host.
import { build } from "esbuild";
import { readFile, writeFile } from "node:fs/promises";

const js = (await build({
  entryPoints: ["src/app.js"], bundle: true, minify: true, format: "iife", target: "es2020", write: false,
})).outputFiles[0].text;
const css = await readFile("app.css", "utf8");
let html = await readFile("index.html", "utf8");

html = html
  .replace('<link rel="stylesheet" href="./app.css">', () => `<style>\n${css}</style>`)
  .replace('<script type="module" src="./src/app.js"></script>', () => `<script>\n${js.replace(/<\/script/gi, "<\\/script")}</script>`);
if (html.includes("./src/app.js") || html.includes("./app.css")) throw new Error("standalone build: placeholders not replaced");

await writeFile("dist/ascii-fold-standalone.html", html);
console.log(`dist/ascii-fold-standalone.html  ${(html.length / 1024).toFixed(1)}kb`);
