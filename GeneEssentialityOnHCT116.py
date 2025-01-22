# ライブラリのimport
import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt

# 各種ファイルのimport
# https://depmap.org/portal/cell_line/ACH-000971?tab=overview 
# 以上DepmapよりHCT116の情報を取得
gene_effect_df = pd.read_csv("gene_effects_ACH-000971.csv") # HCT116で登録されるgene_effect (17334, 5)
gene_effect_df = gene_effect_df.rename(columns={'gene': 'genes'})

# HCT116: https://v3.ogee.info/?#/cellline/large%20intestine/HCT116/summary
# OGEEでHCT116（Source: sanger, Platform: crispr）で登録されている必須遺伝子 
#　tested genes =17995
hct116e_df = pd.read_csv("HCT116_crispr_essential.csv")       # hct116でessentialと登録されているもの (996, 10)

# https://v3.ogee.info/#/downloads
# sanger dataset:	human cell lines and corresponding log fc
sanger_df = pd.read_csv("sanger.txt", sep='\t')               # 全細胞腫で登録されている遺伝子 (5848375, 7)
hct116gene_df = sanger_df[sanger_df["cell_line"]=="HCT116"]   # HCT116で登録されている遺伝子 (17995, 7)

# 結合してEssential/Non-Essentialを付ける
hct116gene_df['essentiality'] = hct116gene_df['genes'].apply(lambda x: 'Essential' if x in hct116e_df['genes'].values else 'Non-Essential')
hct116_pluseffect_df = pd.merge(hct116gene_df, gene_effect_df, on="genes", how="outer")

# 見たい遺伝子がEssentialかどうか
def specific_gene_classify(genes_list):
    temp_df = hct116gene_df[hct116gene_df["genes"].isin(genes_list)]
    return temp_df

def specific_gene_classify_pluseffect(genes_list):
    temp_df = hct116_pluseffect_df[hct116_pluseffect_df["genes"].isin(genes_list)]
    return temp_df

# OGEEに登録されていない遺伝子を除いた。
hct116_pluseffect_df_sort = hct116_pluseffect_df.dropna(subset=['essentiality'])

# show histgram
plt.figure(figsize=(10, 6))
sns.histplot(data=hct116_pluseffect_df_sort, x='gene_effect', hue='essentiality',
             multiple='stack', kde=False, bins=50)

# X軸の範囲を設定
plt.xlim([-2, 1])

# グラフの調整
plt.title('Distribution of Gene Effect by Essentiality')
plt.xlabel('Gene Effect')
plt.ylabel('Count')
plt.grid(True)

plt.show()
