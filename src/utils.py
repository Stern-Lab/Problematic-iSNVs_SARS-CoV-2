import pandas as pd


def assign_transition_type(t):
    transitions = ['AG', 'GA', 'TC', 'CT']
    oxidation = ['CA', 'GT']
    transversions = ['AC', 'TG', 'TA', 'AT', 'GC', 'CG']
    if not isinstance(t, str):
        return 'err'
    if t in transitions:
        return 'ts'
    elif t in oxidation:
        return 'ox'
    elif t[0] == '-':
        return 'ins'
    elif t[1] == '-':
        return 'del'
    elif t in transversions:
        return 'tv'
    else:
        return 'ref'


def add_data_columns(df):
    df['sample'] = df.source.map(lambda s: s.split('-')[0])     # get the sample id
    df['replicate'] = df.source.map(lambda s: s.split('-')[1])  # get replicate id
    df['mutation'] = df['REF'] + df['POS'].astype(int).map(str) + df['ALT']
    df['transition'] = df['REF'] + df['ALT']
    df['type'] = df.transition.map(assign_transition_type)


def merge_variants_w_info(var_df, sym_onset_path, ct_values_path):
    # read info tables
    sym_onset = pd.read_csv(sym_onset_path)
    ct_values = pd.read_csv(ct_values_path)

    # merge
    merged_w_sym_onset = pd.merge(var_df, sym_onset, how='left', on=['sample'])
    merged_variants = pd.merge(merged_w_sym_onset, ct_values[["sample", "N_gene_ct"]], how='left', on=['sample'])

    # drop unnecessary columns & return df
    merged_variants.dropna(how='all', axis='columns')
    return merged_variants
