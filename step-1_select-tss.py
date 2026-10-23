# Step 1b: pick one TSS per gene from the candidate transcripts of
# step-1_annotation-to-transcripts.R, using the ATAC fragments
# (as chd.data.regions.select_tss_from_fragments: the candidate with the most
# fragments within -100..+100 bp of its TSS; ties go to the longest transcript),
# and write transcripts.csv for step 2.
# Usage: python step-1_select-tss.py <candidates.csv> <fragments.tsv.gz> <out transcripts.csv>
# The fragments file must be bgzipped and tabix-indexed (.tbi).
import sys

import pandas as pd
import pysam
import tqdm

if len(sys.argv) != 4:
    sys.exit("usage: python step-1_select-tss.py <candidates.csv> <fragments.tsv.gz> <out transcripts.csv>")
candidates_file, fragments_file, transcripts_file = sys.argv[1:]

candidates = pd.read_csv(candidates_file)

# scaffold genes: dropped by step 2 anyway
scaffold = candidates["chrom"].str.contains("scaffold")
if scaffold.any():
    print(f"dropping {candidates.loc[scaffold, 'ensembl_gene_id'].nunique()} genes on scaffold contigs")
    candidates = candidates.loc[~scaffold].reset_index(drop=True)

# pysam can't fetch from contigs missing in the fragments file (e.g. chloroplast)
contigs = set(pysam.TabixFile(fragments_file).contigs)
missing = ~candidates["chrom"].isin(contigs)
if missing.any():
    print(f"dropping {candidates.loc[missing, 'ensembl_gene_id'].nunique()} genes on contigs without fragments: "
          f"{', '.join(sorted(candidates.loc[missing, 'chrom'].unique()))}")
    candidates = candidates.loc[~missing].reset_index(drop=True)

# Count fragments within -100..+100 bp of each candidate TSS with the same pysam
# fetch as chd's select_tss_from_fragments (chromatinhd 0.4.3), but with the window
# clipped at the contig start: chd passes a negative start for a TSS < 100 bp in,
# which pysam rejects (chd's own fragment loader for step 2 clips the same way).
# chd also breaks ties arbitrarily (unstable sort), so pick here: most fragments,
# then the longest transcript, then transcript id.
fragments_tabix = pysam.TabixFile(fragments_file)
candidates["n_fragments"] = [
    sum(1 for _ in fragments_tabix.fetch(chrom, max(tss - 100, 0), tss + 100))
    for chrom, tss in tqdm.tqdm(zip(candidates["chrom"], candidates["tss"]), total=candidates.shape[0])
]
candidates["length"] = candidates["end"] - candidates["start"] + 1
transcripts = (candidates
    .sort_values(["n_fragments", "length", "transcript"], ascending=[False, False, True], kind="stable")
    .groupby("ensembl_gene_id", sort=False).first()
)
transcripts["ensembl_gene_id"] = transcripts.index

# gene names as step 2 makes them: Glyma.09G186700.Wm82.a4.v1 --> GLYMA_09G186700
transcripts["symbol"] = [g.replace(".Wm82.a4.v1", "").upper().replace(".", "_") for g in transcripts["ensembl_gene_id"]]
transcripts["ensembl_gene_id"] = transcripts["symbol"]
transcripts["transcript_biotype"] = transcripts["gene_biotype"]

columns = ["chrom", "start", "end", "transcript", "strand", "ensembl_gene_id", "tss", "symbol", "transcript_biotype", "n_fragments"]
transcripts[columns].to_csv(transcripts_file, index=False)

print(f"{transcripts.shape[0]} genes; "
      f"{(candidates.groupby('ensembl_gene_id')['tss'].nunique() > 1).sum()} had more than one candidate TSS; "
      f"{(transcripts['n_fragments'] == 0).sum()} have no fragments at the selected TSS")
