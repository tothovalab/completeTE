#!/usr/bin/env python3

import argparse
import pandas as pd

def load_ucsc_rmsk(rmsk_path):
    """
    Load UCSC RepeatMasker annotation
    and return a normalized TE BED-like DataFrame
    """
    colnames = [
        "bin", "swScore", "milliDiv", "milliDel", "milliIns",
        "chr", "start", "end", "genoLeft", "strand",
        "subfamily", "family", "class",
        "repStart", "repEnd", "repLeft", "id"
    ]

    rmsk = pd.read_csv(
        rmsk_path,
        sep="\t",
        header=None,
        names=colnames
    )

    te_bed = rmsk[["chr", "start", "end", "subfamily", "class", "strand"]].copy()
    te_bed = te_bed[~te_bed["subfamily"].str.contains("Simple_repeat|Low_complexity|Satellite", regex=True)]
    te_bed = te_bed[te_bed["strand"].isin(["+", "-"])]

    return te_bed


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Curate complete proviral retrotransposon elements "
            "(LTR–INT–LTR) from UCSC RepeatMasker annotations."
        )
    )

    # other params: tolerance, LTR of interest, int of interest, output folder
    
    parser.add_argument(
        "-r", "--rmsk",
        required=True,
        help="UCSC RepeatMasker annotation file"
    )
    
    parser.add_argument(
        "-i", "--int",
        required=True,
        help="Int subfamily name"
    )
    
    parser.add_argument(
        "-l", "--ltr",
        required=True,
        help="LTR subfamily name"
    )
    
    parser.add_argument(
        "-t", "--tolerance",
        required=False,
        help="tolerance for complete element flanking region",
        default=10000
    )


    return parser.parse_args()

def sort_chromosomes(df):
    # Extract numeric or special parts of chromosome names
    def chr_key(c):
        match = re.match(r"chr(\d+)", c)
        if match:
            return (0, int(match.group(1)))  # autosomes first
        elif c == "chrX":
            return (1, 23)
        elif c == "chrY":
            return (1, 24)
        elif c == "chrM" or c == "chrMT":
            return (2, 25)
        else:
            return (3, c)  # unplaced contigs etc.

    df = df.copy()
    df["chr_sort"] = df["chr"].apply(chr_key)
    df = df.sort_values(["chr_sort", "start"]).drop(columns=["chr_sort"]).reset_index(drop=True)
    return df


def main():
    args = parse_args()

    te_bed = load_ucsc_rmsk(args.rmsk)
    int_of_interest = args.int
    ltr_of_interest = args.ltr
    
    #create log file for these print statements in output dir

    print(f"Loaded {len(te_bed):,} RepeatMasker records")
    print(te_bed.head())
    
    int_elements = te_bed[te_bed['subfamily'] == int_of_interest].sort_values(["chr", "start"]).reset_index(drop=True)
    ltr_elements = te_bed[te_bed['subfamily'] == ltr_of_interest].sort_values(["chr", "start"]).reset_index(drop=True)
    
    int_elements = sort_chromosomes(int_elements)
    ltr_elements = sort_chromosomes(ltr_elements)
    
    tolerance = args.tolerance

    int_elements = int_elements.sort_values(["chr", "start"])
    ltr_elements = ltr_elements.sort_values(["chr", "start"])

    triplets = []

    for _, row in int_elements.iterrows():
        chrom = row["chr"]
        strand = row["strand"]
        int_start = row["start"]
        int_end = row["end"]

        ltrs_chr = ltr_elements[(ltr_elements["chr"] == chrom) & (ltr_elements["strand"] == strand)]
        ltr_up = ltrs_chr[(abs(ltrs_chr["end"] - int_start) <= tolerance) & (ltrs_chr["end"] <= int_start)]
        ltr_down = ltrs_chr[(abs(ltrs_chr["start"] - int_end) <= tolerance) & (ltrs_chr["start"] >= int_end)]

        if not ltr_up.empty and not ltr_down.empty:
            up = ltr_up.iloc[ltr_up["end"].sub(int_start).abs().argsort().iloc[0]]
            down = ltr_down.iloc[ltr_down["start"].sub(int_end).abs().argsort().iloc[0]]

            triplets.append({
                "chr": chrom,
                "strand": strand,
                "ltr_up_start": up["start"],
                "ltr_up_end": up["end"],
                "int_start": int_start,
                "int_end": int_end,
                "ltr_down_start": down["start"],
                "ltr_down_end": down["end"]
            })

    ltr_elements_full = pd.DataFrame(triplets)
    ltr_elements_full

    group_cols = [
        "chr", "strand",
        "ltr_up_start", "ltr_up_end",
        "ltr_down_start", "ltr_down_end"
    ]

    ltr_elements_full_collapsed = (
        ltr_elements_full
        .groupby(group_cols, as_index=False)
        .agg(
            int_start=("int_start", "min"),
            int_end=("int_end", "max")
        )
        .sort_values(["chr", "ltr_up_start", "int_start"])
        .reset_index(drop=True)
    )

    ltr_elements_full_collapsed = ltr_elements_full_collapsed.assign(
        ltr_5_start=np.where(
            ltr_elements_full_collapsed["strand"] == "+",
            ltr_elements_full_collapsed["ltr_up_start"],
            ltr_elements_full_collapsed["ltr_down_start"]
        ),
        ltr_5_end=np.where(
            ltr_elements_full_collapsed["strand"] == "+",
            ltr_elements_full_collapsed["ltr_up_end"],
            ltr_elements_full_collapsed["ltr_down_end"]
        ),
        ltr_3_start=np.where(
            ltr_elements_full_collapsed["strand"] == "+",
            ltr_elements_full_collapsed["ltr_down_start"],
            ltr_elements_full_collapsed["ltr_up_start"]
        ),
        ltr_3_end=np.where(
            ltr_elements_full_collapsed["strand"] == "+",
            ltr_elements_full_collapsed["ltr_down_end"],
            ltr_elements_full_collapsed["ltr_up_end"]
        )
    )

    ltr_elements_full_collapsed.to_csv(out_file)


if __name__ == "__main__":
    main()

