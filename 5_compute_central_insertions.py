import pandas as pd
genomedat = []

pathprefix = "processed/"
from Bio import SeqIO,SeqRecord, Align
for s in SeqIO.parse("reference/XAUV_V0C4.gb","gb"):
    for feat in s.features:
        if feat.type == "CDS":
            genomedat.append({"locus_tag":feat.qualifiers.get("locus_tag","")[0],
                              "cds_start":feat.location.start,
                              "cds_end":feat.location.end,
                              "strand":feat.location.strand})
            locus_tag = feat.qualifiers.get("locus_tag","")[0]
genomedf = pd.DataFrame(genomedat)


def analyze_reads(f):
    strain = f.split("/")[-1].split("_")[0]
    print(strain)
    if "VSF" in f:
        sep = "\t"
    else:
        sep = "\t"
    df = pd.read_csv(f, sep=sep, dtype={"locus_tag":str, 
                                         "cds_start":str,
	                                 "position_mapped":int,
                                         "cds_end":str,
                                         "cds_strand":float,
                                         "cds_product":str}).reset_index()
    df["cds_start"] = [float(str(v).replace("<","")) for v in df.cds_start]
    df["cds_end"] = [float(str(v).replace(">","")) for v in df.cds_end]
    total_insertions = df.shape[0]
    print(f"\tAll sequenced insertions={df.shape[0]}")
    df = df[~(df.locus_tag == "")]
    print(f"\tInsertions in CDS={df.shape[0]}")
    ## Retain only unique insertions
    df = df[["locus_tag","position_mapped"]]\
        .drop_duplicates()\
        .reset_index(drop=True)

    ## Separately, compute central insertions
    dfinserts = df.copy()
    dfinserts = genomedf.merge(dfinserts, on="locus_tag",how="left")
    dfinserts["relative_position"] = dfinserts.position_mapped - dfinserts.cds_start
    dfinserts["cds_size"] = dfinserts.cds_end - dfinserts.cds_start    
    dfinserts = dfinserts[(dfinserts.relative_position > 0.1*dfinserts.cds_size)\
                          & (dfinserts.relative_position < 0.9*(dfinserts.cds_size))]
    dfinserts = dfinserts[["locus_tag","position_mapped"]]\
        .drop_duplicates()\
        .reset_index(drop=True)

    ## compute stats per locus_tag
    df = df.groupby("locus_tag")\
                          .count()\
                          .reset_index()\
                          .rename({"position_mapped":"unique_hit_count"},
                                  axis=1)
    dfinserts = dfinserts.groupby("locus_tag")\
                          .count()\
                          .reset_index()\
                          .rename({"position_mapped":"unique_central_hit_count"},
                                  axis=1)
    mergecols = ["locus_tag"]

    df = genomedf.merge(df, on=mergecols,
                        how="left")
    dfinserts = genomedf.merge(dfinserts, on=mergecols,
                               how="left")
    df.loc[df.unique_hit_count.isna(),"unique_hit_count"] = 0
    df["unique_hit_count"] =df["unique_hit_count"]
    dfinserts.loc[dfinserts.unique_central_hit_count.isna(),"unique_central_hit_count"] = 0
    dfinserts["unique_central_hit_count"] =dfinserts["unique_central_hit_count"] 
    print(f"\tUnique insertions in CDS={df.unique_hit_count.sum()}")
    print(f"\tUnique central insertions in CDS={dfinserts.unique_central_hit_count.sum()}")    

    dffinal = df.merge(dfinserts, on =["locus_tag","cds_start","cds_end","strand"])
    return(dffinal)


flist = [f"{pathprefix}/V0C4_map_all.csv",
         f"{pathprefix}/VSF_map_all.csv"]

v0c4 = analyze_reads(flist[0]).assign(strain = "V0C4")
vsf = analyze_reads(flist[1]).assign(strain = "VSF")
df_full = pd.concat([v0c4, vsf]).reset_index(drop=True)
df_full["locus_id"] = [int(v.split("_")[1]) for v in df_full.locus_tag.values]
df_full.to_csv(f"{pathprefix}/full_vsf_v0c4_read_set.csv",index=False)
