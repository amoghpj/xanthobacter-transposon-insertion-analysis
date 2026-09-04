import pandas as pd

outdir = "processed/"
libraries = ["V0C4_S2_","VSF_S3_"]
suff = "_processed"
for lib in libraries:
    for read in ["R1"]:
        with open(f"{outdir}/{lib}{read}{suff}.csv", "r") as infile:
            with open(f"{outdir}/{lib}{read}{suff}.fasta", "w") as outfile:
                for line in infile:
                    read_id,umi,barcode,genomic_region,genomic_phredscore,prunedGs = line.strip().split(",")
                    outfile.write(f">{read_id.replace('@','')}\n{genomic_region}\n")
