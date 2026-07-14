# biospectral-cli

[![CI](https://github.com/kostyk348/biospectral-cli/actions/workflows/ci.yml/badge.svg)](https://github.com/kostyk348/biospectral-cli/actions/workflows/ci.yml)


Command-line interface for [biospectral](https://github.com/kostyk348/biospectral) —
Dirichlet L-function spectral analysis of DNA and protein sequences.

## Install

```bash
pip install biospectral-cli
```

## Commands

### fingerprint — print the L-spectrum of a sequence

```bash
biospectral-cli fingerprint ATGGTGCACCTGACTCCTGAGG
biospectral-cli fingerprint MKVLSPADKTNVKAAWGKV --kind protein
```

### screen — rank a query against a database

```bash
# built-in reference gene database

[![CI](https://github.com/kostyk348/biospectral-cli/actions/workflows/ci.yml/badge.svg)](https://github.com/kostyk348/biospectral-cli/actions/workflows/ci.yml)

biospectral-cli screen ATGGTGCACCTGACTCCTGTGG

# custom database (FASTA or JSON {"name": "seq"})

[![CI](https://github.com/kostyk348/biospectral-cli/actions/workflows/ci.yml/badge.svg)](https://github.com/kostyk348/biospectral-cli/actions/workflows/ci.yml)

biospectral-cli screen QUERY.fa --db database.fasta --top 10
```

### scan — sliding-window localisation of a reference inside a long sequence

```bash
biospectral-cli scan long_sequence.fa ATGGTGCACCTGACTCCTGAGG
```

### bench — quick built-in benchmark

```bash
biospectral-cli bench
```

Compares the spectral method against Hamming and Levenshtein on a small
synthetic database and prints top-1 accuracy and per-query time.

## How it works

Each sequence is mapped to a real Dirichlet character indexed by position;
its L-function on the critical line `L(t) = sum chi_n n^{-1/2 - it}` is a
fixed-size spectral fingerprint. Comparison reduces to a cheap integral
between fingerprints, so database screening is fast after a one-time
precomputation. See the `biospectral` repository for the full method and
benchmarks.

## License

Apache 2.0
