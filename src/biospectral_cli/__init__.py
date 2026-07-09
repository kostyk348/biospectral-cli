"""biospectral-cli — command-line interface for biospectral.

Commands
--------
fingerprint   Compute and print the L-function spectrum of a sequence.
screen        Rank a query against a sequence database (FASTA / JSON).
scan          Sliding-window localisation of a reference inside a long sequence.
bench         Run a quick built-in benchmark vs Hamming / Levenshtein.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from typing import Dict, List, Tuple

import numpy as np

from biospectral.dna import (
    dna_to_character_seq,
    l_function,
    database_screen,
    sliding_window_scan,
    REFERENCE_GENE_DATABASE,
)
from biospectral.protein import (
    protein_to_character_seq,
    protein_l_function,
    protein_database_screen,
    REFERENCE_PROTEIN_DATABASE,
)


# ---------------------------------------------------------------------------
# FASTA / JSON helpers
# ---------------------------------------------------------------------------

def read_fasta(path: str) -> Dict[str, str]:
    """Minimal FASTA parser -> {header: sequence}."""
    records: Dict[str, str] = {}
    name = None
    buf: List[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
    if name is not None:
        records[name] = "".join(buf)
    return records


def read_db(path: str) -> Dict[str, str]:
    """Read a database from FASTA or JSON ({name: seq})."""
    if path.lower().endswith(".json"):
        with open(path) as fh:
            return json.load(fh)
    return read_fasta(path)


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def cmd_fingerprint(args: argparse.Namespace) -> int:
    seq = args.sequence.upper()
    if args.kind == "dna":
        chi = dna_to_character_seq(seq)
        L = l_function(np.linspace(10.0, 80.0, args.points), chi)
    else:
        chi = protein_to_character_seq(seq)
        L = protein_l_function(np.linspace(10.0, 80.0, args.points), chi)
    np.set_printoptions(precision=4, suppress=True)
    print(f"# {args.kind} sequence ({len(seq)} residues)")
    print(f"# character vector: {chi.tolist()}")
    print(f"# L-function spectrum (Re, Im) at {args.points} points:")
    for re_, im_ in zip(L.real, L.imag):
        print(f"{re_: .6e} {im_: .6e}")
    return 0


def cmd_screen(args: argparse.Namespace) -> int:
    query = args.query.upper()
    if args.db:
        database = read_db(args.db)
    else:
        database = REFERENCE_GENE_DATABASE if args.kind == "dna" else REFERENCE_PROTEIN_DATABASE

    t0 = time.perf_counter()
    if args.kind == "dna":
        results = database_screen(query, database)
    else:
        results = protein_database_screen(query, database)
    dt = (time.perf_counter() - t0) * 1000.0

    print(f"# query ({len(query)} residues) screened against {len(database)} sequences")
    print(f"# time: {dt:.2f} ms\n")
    print(f"{'rank':>4}  {'name':28s} {'residual':>14s}")
    print("-" * 50)
    for i, (name, resid) in enumerate(results[: args.top], 1):
        print(f"{i:>4}  {name:28s} {resid:14.4e}")
    return 0


def cmd_scan(args: argparse.Namespace) -> int:
    long_seq = args.long.upper()
    ref = args.reference.upper()
    t0 = time.perf_counter()
    residuals = sliding_window_scan(long_seq, ref, window_size=len(ref))
    dt = (time.perf_counter() - t0) * 1000.0

    order = np.argsort(residuals)
    print(f"# sliding window: {len(long_seq)} bp, reference {len(ref)} bp")
    print(f"# time: {dt:.2f} ms, {len(residuals)} windows\n")
    print(f"{'pos':>6}  {'residual':>14s}")
    print("-" * 24)
    for pos in order[: args.top]:
        print(f"{pos:>6}  {residuals[pos]:14.4e}")
    return 0


def cmd_bench(args: argparse.Namespace) -> int:
    import random

    def hamming(a, b):
        n = min(len(a), len(b))
        return sum(1 for i in range(n) if a[i] != b[i]) + abs(len(a) - len(b))

    def levenshtein(a, b):
        m, n = len(a), len(b)
        dp = list(range(n + 1))
        for i in range(1, m + 1):
            prev = dp[0]
            dp[0] = i
            for j in range(1, n + 1):
                cur = dp[j]
                dp[j] = min(dp[j] + 1, dp[j - 1] + 1, prev + (a[i - 1] != b[j - 1]))
                prev = cur
        return dp[n]

    rng = random.Random(0)
    bases = "ACGT"
    db = {f"seq_{i}": "".join(rng.choices(bases, k=200)) for i in range(20)}
    queries = [db[f"seq_{i}"][:100] + rng.choice(bases) + db[f"seq_{i}"][101:] for i in range(10)]

    t = np.linspace(10.0, 80.0, 200)

    def spectral(q):
        return database_screen(q, db)

    for name, fn in [("spectral", spectral),
                     ("hamming", lambda q: sorted(((n, hamming(q, s)) for n, s in db.items()), key=lambda r: r[1])),
                     ("levenshtein", lambda q: sorted(((n, levenshtein(q, s)) for n, s in db.items()), key=lambda r: r[1]))]:
        times = []
        correct = 0
        for q, true in zip(queries, [f"seq_{i}" for i in range(10)]):
            t0 = time.perf_counter()
            ranked = fn(q)
            times.append((time.perf_counter() - t0) * 1000.0)
            if ranked[0][0] == true:
                correct += 1
        print(f"{name:12s} top-1 {correct/len(queries)*100:5.1f}%  {np.mean(times):8.2f} ms/query")
    return 0


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="biospectral-cli",
        description="Spectral DNA/protein analysis via Dirichlet L-functions.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    pf = sub.add_parser("fingerprint", help="print the L-spectrum of a sequence")
    pf.add_argument("sequence")
    pf.add_argument("--kind", choices=["dna", "protein"], default="dna")
    pf.add_argument("--points", type=int, default=200)
    pf.set_defaults(func=cmd_fingerprint)

    sc = sub.add_parser("screen", help="rank a query against a database")
    sc.add_argument("query")
    sc.add_argument("--db", help="FASTA or JSON database")
    sc.add_argument("--kind", choices=["dna", "protein"], default="dna")
    sc.add_argument("--top", type=int, default=5)
    sc.set_defaults(func=cmd_screen)

    sw = sub.add_parser("scan", help="sliding-window localisation")
    sw.add_argument("long")
    sw.add_argument("reference")
    sw.add_argument("--top", type=int, default=5)
    sw.set_defaults(func=cmd_scan)

    bn = sub.add_parser("bench", help="quick built-in benchmark")
    bn.set_defaults(func=cmd_bench)

    return p


def main(argv: List[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
