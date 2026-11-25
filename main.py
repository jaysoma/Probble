#!/usr/bin/env python3
"""
Probble (Protein + Wobble) — AlphaFold structure library + protein renderer.

Landing (catalog.html) is the browse-by-fold-class library; every card opens explore.html,
the four-tab renderer (Explore / Assemble / Fold / Atoms). It runs on a locally-built
3Dmol.js (3Dmol-dev.js) carrying the per-atom sphere opacity feature (upstream issue #166):
the "spheres · confidence" representation maps each residue's pLDDT (the B-factor column of
every AlphaFold model) onto per-atom opacity. The DPAM "step the paper" walkthrough this app
began as was stripped 2026-08-09 — git 41caa40 has it whole. viewer.html and alphatest.html
are unlinked dev/test pages.

A static app behind the portfolio gateway (GithubRoot/app.py). The gateway
reverse-proxies /foldview/* -> 127.0.0.1:<PORT>/* (prefix stripped), so every
page is served from this folder's root and all in-page paths stay relative.

Port comes from the launcher via $env:PORT (start-all.ps1); falls back to 8005.
Landing -> catalog.html (the browse-by-fold-class library).
"""
import http.server
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
os.chdir(ROOT)
PORT = int(os.environ.get("PORT", 8005))
HOST = "127.0.0.1"


class Handler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header("Cache-Control", "no-store, must-revalidate")
        self.send_header("Access-Control-Allow-Origin", "*")
        super().end_headers()

    def do_GET(self):
        # gateway hits "/" (route root) -> the library is the front door
        if self.path in ("", "/"):
            self.path = "/catalog.html"
        return super().do_GET()

    def log_message(self, fmt, *args):
        try:
            print(f"  {args[0]}  {args[1]}", flush=True)
        except Exception:
            pass


if __name__ == "__main__":
    print(f"\n  FoldView  ->  http://{HOST}:{PORT}   (landing: catalog.html)\n", flush=True)
    http.server.ThreadingHTTPServer((HOST, PORT), Handler).serve_forever()
