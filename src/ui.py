"""
`maillard ui`: one page on this machine that takes a spec and returns the HTML report.

Standard library only (http.server). The page posts the spec text and the verb to /run; the server
validates it, runs the same functions the verbs use (src.api), and returns the report that
`--report` would have written. Nothing leaves the machine; there is no account and no state.
"""
from __future__ import annotations

import json
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any, Dict, Optional, Tuple

import yaml

PAGE = """<!doctype html>
<meta charset="utf-8"><title>Maillard</title>
<style>
 body{font:15px/1.45 system-ui,sans-serif;max-width:60rem;margin:2rem auto;padding:0 1rem;color:#1e2a2c}
 textarea{width:100%;height:24rem;font:13px/1.4 ui-monospace,monospace;padding:.6rem;border:1px solid #b9c4c1;border-radius:6px}
 input,select,button{font:inherit;padding:.4rem .7rem;border:1px solid #b9c4c1;border-radius:6px;background:#fff}
 button{background:#2b5da8;color:#fff;border:none;cursor:pointer} button:disabled{opacity:.5}
 .row{display:flex;gap:.8rem;align-items:center;flex-wrap:wrap;margin:.8rem 0}
 #out{margin-top:1.2rem;border-top:1px solid #d6dbd9;padding-top:1rem} .err{color:#b23a3a;white-space:pre-wrap}
 small{color:#5e6b6e}
</style>
<h1>Maillard</h1>
<p>Paste a spec (the same YAML the command line takes), pick the verb, run. The result is the report
<code>--report</code> writes. Absolute numbers are unreliable; read the ratios and the refusals.</p>
<textarea id="spec">__TEMPLATE__</textarea>
<div class="row">
 <select id="verb"><option value="compare">compare (two arms, a and b)</option><option value="predict">predict (arm a)</option></select>
 <label>calibration file <input id="cal" placeholder="results/user/my_lab/calibration_....json (optional)" size="46"></label>
 <button id="run">Run</button> <small id="status"></small>
</div>
<div id="out"></div>
<script>
const run=document.getElementById('run'), out=document.getElementById('out'), st=document.getElementById('status');
run.onclick=async()=>{run.disabled=true; st.textContent='running...'; out.innerHTML='';
 try{const r=await fetch('/run',{method:'POST',headers:{'content-type':'application/json'},
   body:JSON.stringify({verb:document.getElementById('verb').value, spec:document.getElementById('spec').value, calibration:document.getElementById('cal').value||null})});
  const t=await r.text(); if(!r.ok){out.innerHTML='<div class="err">'+t.replace(/</g,'&lt;')+'</div>'} else {const f=document.createElement('iframe'); f.style.width='100%'; f.style.height='70vh'; f.style.border='1px solid #d6dbd9'; f.srcdoc=t; out.appendChild(f)}
 }catch(e){out.innerHTML='<div class="err">'+e+'</div>'} st.textContent=''; run.disabled=false;};
</script>
"""


def run_request(body: Dict[str, Any]) -> Tuple[int, str, str]:
    """(status, content type, text) for a posted {verb, spec, calibration}."""
    from src import api
    from src.comparative_cli import SpecError
    from src.report_html import render_compare_report, render_predict_report

    verb = str(body.get("verb") or "compare")
    try:
        document = yaml.safe_load(str(body.get("spec") or "")) or {}
        if not isinstance(document, dict):
            raise SpecError("the spec must be a YAML mapping")
        calibration = body.get("calibration") or None
        if verb == "compare":
            if not ({"a", "b"} <= set(document)):
                raise SpecError("compare needs two arms, 'a' and 'b'")
            payload = api.compare(document["a"], document["b"], calibration=calibration)
            return 200, "text/html; charset=utf-8", render_compare_report(payload)
        if verb == "predict":
            spec = document.get("a") if "a" in document else document
            payload = api.predict(spec, calibration=calibration)
            return 200, "text/html; charset=utf-8", render_predict_report(payload)
        raise SpecError(f"unknown verb {verb!r}")
    except (SpecError, ValueError, KeyError, yaml.YAMLError, FileNotFoundError) as exc:
        return 400, "text/plain; charset=utf-8", f"{type(exc).__name__}: {exc}"


class _Handler(BaseHTTPRequestHandler):
    def log_message(self, *args):  # quiet
        return

    def _send(self, status: int, ctype: str, text: str) -> None:
        data = text.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        from src.comparative_cli import SPEC_TEMPLATE

        if self.path in ("/", "/index.html"):
            self._send(200, "text/html; charset=utf-8", PAGE.replace("__TEMPLATE__", SPEC_TEMPLATE.replace("<", "&lt;")))
        else:
            self._send(404, "text/plain; charset=utf-8", "not found")

    def do_POST(self):
        if self.path != "/run":
            self._send(404, "text/plain; charset=utf-8", "not found")
            return
        length = int(self.headers.get("Content-Length") or 0)
        try:
            body = json.loads(self.rfile.read(length).decode("utf-8") or "{}")
        except json.JSONDecodeError:
            self._send(400, "text/plain; charset=utf-8", "the request body must be JSON")
            return
        self._send(*run_request(body))


def make_server(port: int = 8765, host: str = "127.0.0.1") -> ThreadingHTTPServer:
    return ThreadingHTTPServer((host, port), _Handler)


def serve(port: int = 8765, open_browser: bool = True) -> int:
    server = make_server(port)
    url = f"http://127.0.0.1:{server.server_address[1]}/"
    print(f"maillard ui at {url}  (Ctrl-C to stop)")
    if open_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
    return 0
