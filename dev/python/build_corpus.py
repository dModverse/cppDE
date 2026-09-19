"""Merge the tap log into corpus.json.

Usage: python dev/python/build_corpus.py [logdir] [out]

Output: {"ode", "event", "rootfunc", "fun", "cvode"} -> list of distinct call
kwargs, in first-seen order.
"""

import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    logdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "tap_log")
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "corpus.json")
    seen = set()
    corpus = {"ode": [], "event": [], "rootfunc": [], "fun": [], "cvode": []}
    kinds = {
        ("codegen_cppODE", "generate_ode_cpp"): "ode",
        ("codegen_cppODE", "generate_event_code"): "event",
        ("codegen_cppODE", "generate_rootfunc_code"): "rootfunc",
        ("codegen_cppFUN", "generate_fun_cpp"): "fun",
        ("codegen_cvode", "generate_cvode_cpp"): "cvode",
    }
    drop = {"outdir", "modelname", "version", "srcfile"}
    for path in sorted(glob.glob(os.path.join(logdir, "*.jsonl"))):
        with open(path, encoding="utf-8") as f:
            for line in f:
                try:
                    rec = json.loads(line)
                except ValueError:
                    continue
                kind = kinds.get((rec["module"], rec["func"]))
                if kind is None:
                    continue
                kw = {k: v for k, v in rec["kwargs"].items() if k not in drop}
                key = json.dumps([kind, kw], sort_keys=True)
                if key in seen:
                    continue
                seen.add(key)
                corpus[kind].append(kw)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(corpus, f, indent=0, sort_keys=True)
        f.write("\n")
    print(out, {k: len(v) for k, v in corpus.items()})


if __name__ == "__main__":
    main()
