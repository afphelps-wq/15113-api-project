"""Build the small demo dataset committed to the repo.

Source: GEO series GSE205154 (289 bulk RNA-seq PDAC tumors, FFPE, Homo sapiens).
Download the full files into test_data/ first (they are git-ignored because of size):
  GSE205154_Gene_Level_Counts_Estimates.txt   (supplementary file on the GEO page)
  GSE205154_series_matrix.txt                  (metadata)

Output (demo_data/):
  demo_counts.csv    genes x samples integer counts, plus gene_id and gene_symbol columns
  demo_metadata.csv  sample, tumor_type, tissue

Run from the repo root:  python scripts/make_demo_data.py
"""
from pathlib import Path

import pandas as pd

SEED = 0
N_PER_GROUP = 30          # 30 Primary + 30 Met
MIN_COUNT = 10            # a gene must have >= MIN_COUNT reads in >= MIN_SAMPLES samples
MIN_SAMPLES = 6           # 10% of the 60 demo samples

SRC = Path("test_data")
OUT = Path("demo_data")


def read_series_matrix(path):
    """Return a per-sample metadata table parsed from a GEO series matrix file."""
    titles, chars = None, {}
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        if parts[0] == "!Sample_title":
            titles = [x.strip('"') for x in parts[1:]]
        elif parts[0] == "!Sample_characteristics_ch1":
            values = [x.strip('"') for x in parts[1:]]
            name = values[0].split(":")[0].strip().replace(" ", "_")
            chars[name] = [v.split(":", 1)[1].strip() for v in values]
    return pd.DataFrame(chars, index=titles).rename_axis("sample")


def main():
    meta = read_series_matrix(SRC / "GSE205154_series_matrix.txt")
    raw = pd.read_csv(SRC / "GSE205154_Gene_Level_Counts_Estimates.txt", sep="\t", low_memory=False)

    picked = pd.concat([
        meta[meta["tumor_type"] == group].sample(N_PER_GROUP, random_state=SEED)
        for group in ("Primary", "Met")
    ])
    counts = raw[picked.index].round().astype(int)     # GEO gives estimated (non-integer) counts

    keep = (counts >= MIN_COUNT).sum(axis=1) >= MIN_SAMPLES
    # Gene label: HGNC symbol, else GEO gene_name, else the Ensembl ID
    symbol = raw["hgnc_symbol"].fillna(raw["gene_name"]).fillna(raw["ID"])
    out = pd.concat([raw[["ID"]].rename(columns={"ID": "gene_id"}),
                     symbol.rename("gene_symbol"), counts], axis=1)[keep]

    OUT.mkdir(exist_ok=True)
    out.to_csv(OUT / "demo_counts.csv", index=False)
    picked[["tumor_type", "tissue"]].to_csv(OUT / "demo_metadata.csv")
    print(f"demo_counts.csv: {out.shape[0]} genes x {counts.shape[1]} samples")
    print(picked.groupby(["tumor_type", "tissue"]).size())


if __name__ == "__main__":
    main()
