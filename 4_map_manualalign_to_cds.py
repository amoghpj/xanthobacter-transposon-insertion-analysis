"""
conda env: seq 
"""
from Bio import SeqIO
import pandas as pd
from intervaltree import IntervalTree
refpath = "reference/XAUV_V0C4.gb"
pathprefix= "processed"
mapped_paths = [f"{pathprefix}/V0C4_S2_R1_processed_mapped.sam",
                f"{pathprefix}/VSF_S3_R1_processed_mapped.sam"]


Features = []
tree = IntervalTree()
for s in SeqIO.parse(refpath, "gb"):
    for feature in s.features:
        if feature.type == "CDS":
            Features.append(feature)
            tree[feature.location.start:feature.location.end] = feature


def get_write_location(pos, fname, tree):
    hits = tree[pos]
    if len(hits) > 0:
        for _fdat in hits:
            _f = _fdat.data
            f = _f.qualifiers
            with open(outname,"a+") as outfile:
                outfile.write(f"{f.get('locus_tag')[0]}\t{_f.location.start}\t{pos}\t{_f.location.end}\t{_f.location.strand}\t{f['product'][0]}\n")                
for mpath in mapped_paths:
    fname = mpath.split("/")[-1]
    print(fname)
    outname = f"{pathprefix}/{fname}_map_details.csv"
    with open(outname,"w") as outfile:
        outfile.write("locus_tag\tcds_start\tposition_mapped\tcds_end\tcds_strand\tcds_product\n")
    with open(mpath ,"r") as infile:
        for line in infile:
            isheader = False
            for s in ["@HD","@SQ","@PG"]: ## headers
                if (s in line):
                    isheader = True
            if not isheader:
                align = line.strip().split("\t")
                if int(align[1]) != 4:
                    get_write_location(int(align[3]), fname, tree)
    pd.read_csv(outname,sep="\t")\
      .sort_values(by="locus_tag")\
      .to_csv(outname, sep="\t",index=False)
                    
