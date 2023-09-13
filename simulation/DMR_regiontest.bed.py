import pandas as pd

df = pd.read_table('../metilene_DMR_simulation_scripts/bg/DMR_annotation_beta_40_3_1.2_0.8_ratio.0.6.bed', header=None).sort_values(1)

df[3] = [0]+list(df[2][:-1])

df = df.loc[(df[1]>df[3])|(df[2]>df[3])]

df[[0,1,2]].to_csv('../metilene_DMR_simulation_scripts/bg/DMR_annotation_beta_40_3_1.2_0.8_ratio.0.6.short.bed', sep='\t',index=False,header=False)
