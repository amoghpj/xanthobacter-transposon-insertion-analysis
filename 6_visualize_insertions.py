import os
import sys
from tqdm import tqdm
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Bio import SeqIO,SeqRecord, Align
from scipy.stats import norm
from adjustText import adjust_text
from matplotlib.ticker import MaxNLocator
from statistics import NormalDist

if not os.path.exists("./img"):
    os.mkdir("./img/")
pathprefix = "processed"
sns.set(
        font_scale=1.5,
        style = "ticks")
### Step 0 -- read the genome annotations --- TODO AMOGH - change  here with PGAP annotation
genomedat = []
for s in SeqIO.parse("reference/XAUV_V0C4.gb","gb"):
    for feat in s.features:
        if feat.type == "CDS":
            genomedat\
                .append({"locus_tag":feat.qualifiers.get("locus_tag","")[0],
                              "cds_start":feat.location.start,
                              "cds_end":feat.location.end,
                         "product":feat.qualifiers.get("product",[])[0],
                              "strand":feat.location.strand})
            locus_tag = feat.qualifiers.get("locus_tag","")[0]
genomedf = pd.DataFrame(genomedat)
genomedf["locus_id"] = [int(v.split("_")[1]) for v in genomedf.locus_tag.values]
####

#### Step 1.1 Use full data
fulldf = pd.read_csv(f"{pathprefix}/full_vsf_v0c4_read_set.csv")

centhits = fulldf.pivot(index="locus_id",columns="strain",
                    values="unique_central_hit_count").reset_index()
centhits["|VSF|/|V0C4|"] = (centhits.VSF)/(centhits.V0C4) # READ_CORRECTION_FACTOR* 
centhits = centhits.merge(genomedf, on="locus_id")
df = centhits[["locus_id","|VSF|/|V0C4|","VSF","V0C4","product"]].set_index("locus_id")

## Drop all genes where number of insertions is less than 4 for each.
INSERTIONS_CUTOFF = 8
print("pre insertion filter", df.shape)
# suffix = "triangle"
# condition = (df.V0C4 + df.VSF ) < INSERTIONS_CUTOFF
suffix = "square"
INS_CUTOFF_VSF = 10
INS_CUTOFF_V0C4 = 8
condition = (df.V0C4 < INS_CUTOFF_V0C4) & (df.VSF  < INS_CUTOFF_VSF)
dflowcounts = df[condition]
df = df[~condition]
print("post insertion filter", df.shape)

READ_CORRECTION_FACTOR = df.V0C4.sum()/df.VSF.sum()

df["|VSF|/|V0C4|"] = READ_CORRECTION_FACTOR* (df.VSF)/(df.V0C4)
print("Insertion correction factor",READ_CORRECTION_FACTOR)

###############
### Step 2 -- Merge counts with annotations
df["Locus tag"] = "LLIMPF_" + df.index.astype(str).str.zfill(5)
df["log2(|VSF|/|V0C4|)"] =np.log2(df["|VSF|/|V0C4|"])
dfasym = df[np.isinf(df["log2(|VSF|/|V0C4|)"].values)]
df = df[~np.isinf(df["log2(|VSF|/|V0C4|)"].values)]

median = np.median(df["|VSF|/|V0C4|"])
std_dev = np.std(df["|VSF|/|V0C4|"])
df["Log2(VSF/V0C4)"] = np.log2(df["|VSF|/|V0C4|"])

### compute MAD
median = np.median(df["log2(|VSF|/|V0C4|)"])
print("median", median)
MAD = np.median(abs(df["log2(|VSF|/|V0C4|)"] - median))
sigma = 1.4826 * MAD
df["z"] = (df["log2(|VSF|/|V0C4|)"] - median)/sigma
df["p"] = 2 * (1 - norm.cdf(abs(df["z"])))

### FDR
from scipy.stats import false_discovery_control
df["q"] = false_discovery_control(df["p"], method="bh")

df.loc[df["p"] == 0, "p"] = np.nan

df["-log10(p)"] = -np.log10(df["p"])
df.loc[df["-log10(p)"].isna(), "-log10(p)"] = 16


### Bonferroni
alpha = 0.05
numgenes = df.shape[0]
print(numgenes)
pthresh = alpha/numgenes

### 2.1 Visualize Bonferroni corrected p-values
fc_pvals = df[["log2(|VSF|/|V0C4|)","-log10(p)"]].sort_values(by="-log10(p)", ascending=False)
fc_cutoff = abs(fc_pvals[fc_pvals["-log10(p)"] > -np.log10(pthresh)]["log2(|VSF|/|V0C4|)"].values[-1])

df = df.assign(is_significant_FWER = 0)
df.loc[abs(df["log2(|VSF|/|V0C4|)"]) > fc_cutoff, "is_significant_FWER"] = 1
high = df[df["|VSF|/|V0C4|"] > 2**fc_cutoff]
low = df[df["|VSF|/|V0C4|"] < 2**(-fc_cutoff)]

g = sns.relplot(data=df,
            x="log2(|VSF|/|V0C4|)",color="k",s=10,alpha=0.2,height=6,aspect=1.5,
            y="-log10(p)",edgecolor=None)
ax = g.axes.flatten()[0]

annotdf = df[(((df["log2(|VSF|/|V0C4|)"]) > fc_cutoff) |\
             ((df["log2(|VSF|/|V0C4|)"]) < -fc_cutoff))\
             & (df["-log10(p)"] > -np.log10(pthresh))]


#annotdf = df[(df["-log10(p)"] > -np.log10(pthresh))]
print("Bonferroni corrected", annotdf.shape)
texts = []
for i, row in annotdf.iterrows():
    ax.plot(row["log2(|VSF|/|V0C4|)"], row["-log10(p)"], 'ro',markeredgecolor="k")
    texts.append(ax.text(row["log2(|VSF|/|V0C4|)"] , row["-log10(p)"], row["Locus tag"], fontsize=8))
p = Rectangle((-fc_cutoff, 0 ),2*fc_cutoff, -np.log10(pthresh),alpha=0.1,color="k")
ax.add_patch(p)
adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))
plt.savefig("img/volcano-plot_full_bonf.pdf")
plt.close("all")



### Step 3 -- the full scatter
fig = plt.figure(figsize=(18,6), constrained_layout=True)
axes = fig.subplot_mosaic("AAAAAB")
sns.set(font_scale=1.5,style="ticks")

g = sns.scatterplot(data=df,
            x="locus_id",color="k",s=10,alpha=0.2,
            y="|VSF|/|V0C4|",edgecolor=None, ax=axes["A"]) # aspect=3

ax = axes["A"] #g.axes.flatten()[0]
for i, row in high.iterrows():
    ax.plot(i, row["|VSF|/|V0C4|"], 'ro',ms=8,mec="k")
for i, row in low.iterrows():
    ax.plot(i, row["|VSF|/|V0C4|"], 'ro',ms=8,mec="k")


for i, row in dfasym.iterrows():
    y = None
    if (row["VSF"] > 0) and (row["V0C4"] == 0):
        y = 256
    if (row["VSF"] == 0) and (row["V0C4"]> 0):
        y = 1./256.
    if y is not None:
        ax.plot(i, y, 'ro',ms=8,mec="k")

color = "#329ddb"
color = "r"
for i, row in high.iterrows():
    if (i >= 24955 ) and (i<= 25115):
        ax.plot(i, row["|VSF|/|V0C4|"], 'o',c=color,ms=8)
for i, row in low.iterrows():
    if (i >= 24955 ) and (i<= 25115):
        ax.plot(i, row["|VSF|/|V0C4|"], 'o',c=color,ms=8)

p1 = Rectangle((24955,1/250.), 250, 250, color="b", alpha=0.1)
ax.add_patch(p1)
ax.axhline(2**fc_cutoff , color='k')
ax.axhline(2**(-fc_cutoff), color='k')
ax.set_yscale("log",base=2)
ax.set_xlim(0,25835)
ax.set_ylim(2.**(-8), 2.**(8))

ax.xaxis.set_major_locator(MaxNLocator(nbins=20))
ax.axhline(1, alpha=0.4, color="k")

g = sns.histplot(data=df, y="Log2(VSF/V0C4)", color="k",
                 bins=50, ax=axes["B"])

ax = axes["B"]
ax.set_ylim(-7,7)
X = np.linspace(-7, 7, 100)

ax.axhline(fc_cutoff, color='k', alpha=0.5)
ax.axhline(-fc_cutoff, color='k', alpha=0.5)



ax.set_yticks([])
ax.set_ylabel("")
plt.tight_layout()
plt.locator_params(axis='x', nbins=5)

plt.savefig(f"img/cent-hits-{suffix}_bonf_perstraincutoff.pdf")
plt.close("all")


### 3.2 Full scatter with FDR
### 2.2 Visualize FDR corrected p-values
fc_qvals = df[["log2(|VSF|/|V0C4|)","q"]].sort_values(by="q", ascending=False)
print("Filtered based on qvals:")
fc_qvalsfilt = fc_qvals[fc_qvals["q"] < 0.05]
fc_cutoff = abs(fc_qvalsfilt["log2(|VSF|/|V0C4|)"].values[0])
print(fc_qvalsfilt.shape, fc_cutoff)
df["-log10(q)"] = -np.log10(df["q"])

high = df[df["|VSF|/|V0C4|"] > 2**fc_cutoff]
low = df[df["|VSF|/|V0C4|"] < 2**(-fc_cutoff)]

g = sns.relplot(data=df,
            x="log2(|VSF|/|V0C4|)",color="k",s=10,alpha=0.2,height=6,aspect=1.5,
            y="-log10(q)",edgecolor=None)
ax = g.axes.flatten()[0]

annotdf_fdr = df[(abs(df["log2(|VSF|/|V0C4|)"]) > fc_cutoff) & (df["q"] < 0.05)]

print("Pass Bonferroni but missing in FDR:", [v for v in list(annotdf.index) if v not in list(annotdf_fdr.index)] )
print("q value cutoff < 0.05", annotdf_fdr.shape)
texts = []
for i, row in annotdf_fdr.iterrows():
    ax.plot(row["log2(|VSF|/|V0C4|)"], row["-log10(q)"], 'ro',markeredgecolor="k")
    texts.append(ax.text(row["log2(|VSF|/|V0C4|)"] , row["-log10(q)"], row["Locus tag"], fontsize=8))
p = Rectangle((-fc_cutoff, 0 ),2*fc_cutoff, -np.log10(0.05),alpha=0.1,color="k")
ax.add_patch(p)
plt.tight_layout()
plt.savefig("img/volcano-plot_full_fdr.pdf")
plt.savefig("img/volcano-plot_full_fdr.png")
plt.close("all")


### Step 3 -- the full scatter
fig = plt.figure(figsize=(18,6), constrained_layout=True)
axes = fig.subplot_mosaic("AAAAAB")
sns.set(font_scale=1.5,style="ticks")

g = sns.scatterplot(data=df,
            x="locus_id",color="k",s=10,alpha=0.2,
            y="|VSF|/|V0C4|",edgecolor=None, ax=axes["A"]) # aspect=3

ax = axes["A"] #g.axes.flatten()[0]
for i, row in high.iterrows():
    ax.plot(i, row["|VSF|/|V0C4|"], 'ro',ms=8,mec="k")
for i, row in low.iterrows():
    ax.plot(i, row["|VSF|/|V0C4|"], 'ro',ms=8,mec="k")


for i, row in dfasym.iterrows():
    y = None
    if (row["VSF"] > 0) and (row["V0C4"] == 0):
        y = 256
    if (row["VSF"] == 0) and (row["V0C4"]> 0):
        y = 1./256.
    if y is not None:
        ax.plot(i, y, 'ro',ms=8,mec="k")

color = "#329ddb"
color = "r"
for i, row in high.iterrows():
    if (i >= 24955 ) and (i<= 25115):
        ax.plot(i, row["|VSF|/|V0C4|"], 'o',c=color,ms=8)
for i, row in low.iterrows():
    if (i >= 24955 ) and (i<= 25115):
        ax.plot(i, row["|VSF|/|V0C4|"], 'o',c=color,ms=8)
# print(fc_cutoff)
p1 = Rectangle((24955,1/250.), 250, 250, color="b", alpha=0.1)
ax.add_patch(p1)
ax.axhline(2**fc_cutoff , color='k')
ax.axhline(2**(-fc_cutoff), color='k')
ax.set_yscale("log",base=2)
ax.set_xlim(0,25835)
ax.set_ylim(2.**(-8), 2.**(8))

ax.xaxis.set_major_locator(MaxNLocator(nbins=20))
ax.axhline(1, alpha=0.4, color="k")

g = sns.histplot(data=df, y="Log2(VSF/V0C4)", color="k",
                 bins=50, ax=axes["B"])

ax = axes["B"]
ax.set_ylim(-8,8)
X = np.linspace(-7, 7, 100)

ax.axhline(fc_cutoff, color='k', alpha=0.5)
ax.axhline(-fc_cutoff, color='k', alpha=0.5)

ax.set_yticks([])
ax.set_ylabel("")
plt.tight_layout()
plt.locator_params(axis='x', nbins=5)

plt.savefig(f"img/cent-hits-{suffix}_fdr_perstraincutoff.pdf")
plt.savefig(f"img/cent-hits-{suffix}_fdr_perstraincutoff.png")
plt.close("all")

df = df.assign(is_significant_FDR = 0)
df.loc[abs(df["log2(|VSF|/|V0C4|)"]) > fc_cutoff, "is_significant_FDR"] = 1

df[["Locus tag","product","Log2(VSF/V0C4)","z","-log10(p)","-log10(q)","is_significant_FWER","is_significant_FDR"]].to_csv(f"{pathprefix}/processed_insertions.csv",index=False)
