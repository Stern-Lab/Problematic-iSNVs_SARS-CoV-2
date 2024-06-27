import pandas as pd


def merge_variants_w_info(var_df, sym_onset_path="data/sym_onset.csv", ct_values_path="data/ct_values.csv"):
    # read info tables
    sym_onset = pd.read_csv(sym_onset_path)
    ct_values = pd.read_csv(ct_values_path)

    # merge
    merged_w_sym_onset = pd.merge(var_df, sym_onset, how='left', on=['sample'])
    merged_variants = pd.merge(merged_w_sym_onset, ct_values[["sample", "N_gene_ct"]], how='left', on=['sample'])

    # drop unnecessary columns & return df
    merged_variants.dropna(how='all', axis='columns')
    return merged_variants