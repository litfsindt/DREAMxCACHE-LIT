import pandas as pd

#read parquet
df= pd.read_parquet(f"/projects/cache_challenge/DREAM_CACHE/datasets/PGK2_selection/Deduplicated_data/PGK2_selection_deduplicated.parquet", engine="fastparquet")

#filter using given rules
filtered_df = df[(df["count_PGK2"] >= 3) & (df["count_NTC"] ==0) & (df["historic_hits"] < 5) & (df["count_PGK2_with_inhibitor"] == 0)].reset_index(drop=True)

filtered_df[["SMILES", "compound"]].to_csv("actives.smi",sep="\t",index=False,header=False)
