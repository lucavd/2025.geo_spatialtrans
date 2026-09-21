"""tools/R2b_run_vlm.py — R2b: invia i 24 ritagli a un modello visione-linguaggio OpenAI-compatibile con il prompt pre-registrato.

Una immagine per richiesta, contesto pulito, temperature 0, una sola richiesta per finestra (nessun retry sul contenuto: se il
modello non risponde con JSON valido, la risposta grezza e' conservata e la finestra risulta non valutata).
La chiave API e' letta dalla variabile d'ambiente indicata (mai stampata, mai salvata). Le risposte grezze e i metadati
(id modello restituito, usage, timestamp) sono salvati per finestra: sono la registrazione esatta della corsa.
Uso: python R2b_run_vlm.py --base-url https://api.deepseek.com --model deepseek-flash --key-env DEEPSEEK_API_KEY \
        --images <cartella con i PNG e manifest.csv> --prompt <R2b_model_prompt.md> --out <cartella output> [--detail high]
Output: <out>/<win_id>.json (JSON del modello, se parsabile), <out>/<win_id>.raw.txt (testo grezzo), <out>/run_log.csv
"""
import sys, os, re, json, base64, time, argparse, csv
from pathlib import Path
import requests

def prompt_block(md_text):
    lines = md_text.splitlines(); seps = [i for i, l in enumerate(lines) if l.strip() == "====="]
    assert len(seps) >= 2, "righe separatrici ===== non trovate"
    return "\n".join(lines[seps[0] + 1:seps[1]]).strip()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-url", required=True); ap.add_argument("--model", required=True); ap.add_argument("--key-env", required=True)
    ap.add_argument("--images", required=True); ap.add_argument("--prompt", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--detail", default=None); ap.add_argument("--only", nargs="*"); ap.add_argument("--max-tokens", type=int, default=8000)
    ap.add_argument("--thinking", choices=["enabled", "disabled", "default"], default="default")
    a = ap.parse_args()
    key = os.environ.get(a.key_env); assert key, f"variabile {a.key_env} assente"
    prompt = prompt_block(Path(a.prompt).read_text())
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    imgs = sorted(Path(a.images).glob("*.png"))
    if a.only: imgs = [p for p in imgs if p.stem in a.only]
    logf = out / "run_log.csv"; new = not logf.exists()
    with open(logf, "a", newline="") as lf:
        w = csv.writer(lf)
        if new: w.writerow(["win_id", "model_requested", "thinking", "max_tokens", "model_returned", "status", "http", "n_nuclei", "n_points", "prompt_tokens", "completion_tokens", "reasoning_tokens", "latency_s", "finish_reason", "timestamp"])
        for p in imgs:
            if (out / f"{p.stem}.json").exists(): print(p.stem, "gia' fatto"); continue
            b64 = base64.b64encode(p.read_bytes()).decode()
            img_block = {"type": "image_url", "image_url": {"url": f"data:image/png;base64,{b64}"}}
            if a.detail: img_block["image_url"]["detail"] = a.detail
            body = {"model": a.model, "temperature": 0, "max_tokens": a.max_tokens,
                    **({"thinking": {"type": a.thinking}} if a.thinking != "default" else {}),
                    "messages": [{"role": "user", "content": [{"type": "text", "text": prompt + f"\n\nThe image file name is: {p.name}"}, img_block]}]}
            t0 = time.time()
            try:
                r = requests.post(a.base_url.rstrip("/") + "/chat/completions", headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json"},
                                  json=body, timeout=600)
                http = r.status_code; data = r.json()
            except Exception as e:
                http = -1; data = {"error": str(e)}
            lat = time.time() - t0
            txt = ""; usage = data.get("usage", {}) if isinstance(data, dict) else {}; model_ret = data.get("model", "") if isinstance(data, dict) else ""
            finish = ""
            try:
                ch = data["choices"][0]; txt = ch["message"]["content"] or ""; finish = ch.get("finish_reason", "")
            except Exception: pass
            (out / f"{p.stem}.raw.txt").write_text(txt if txt else json.dumps({k: v for k, v in data.items() if k != "choices"}, indent=1))
            json.dump(data, open(out / f"{p.stem}.resp.json", "w"), indent=1)          # risposta completa (usage, reasoning, ...)
            status, n, npts = "no_json", "", ""
            m = re.search(r"\{.*\}", txt, re.S)
            if m:
                try:
                    obj = json.loads(m.group(0)); obj.setdefault("image", p.name); obj["image"] = p.name
                    json.dump(obj, open(out / f"{p.stem}.json", "w"), indent=1); status = "ok"
                    n = obj.get("n_nuclei"); npts = len(obj.get("points") or [])
                except json.JSONDecodeError: status = "bad_json"
            w.writerow([p.stem, a.model, a.thinking, a.max_tokens, model_ret, status, http, n, npts, usage.get("prompt_tokens"), usage.get("completion_tokens"), (usage.get("completion_tokens_details") or {}).get("reasoning_tokens"), round(lat, 1), finish, time.strftime("%Y-%m-%dT%H:%M:%S")]); lf.flush()
            print(f"{p.stem} {status} http={http} n={n} pts={npts} tok_in={usage.get('prompt_tokens')} tok_out={usage.get('completion_tokens')} {lat:.0f}s finish={finish}")

if __name__ == "__main__":
    main()
