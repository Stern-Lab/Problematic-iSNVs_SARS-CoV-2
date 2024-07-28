import pandas as pd
import filtering_functions
import utils

# TODO: add each function's description

def get_common_muts(df, replicates=True):
    # group by variant definition and count unique sources
    if replicates:
        isnv_counts = df.groupby(['POS', 'REF', 'ALT', 'mutation'])['sample'].nunique().reset_index(name='source_count')

    else:
        isnv_counts = df.groupby(['POS', 'REF', 'ALT', 'mutation'])['source'].nunique().reset_index(name='source_count')

    return isnv_counts


def find_recurrent_mutations_in_variants(df, variant_sample_counts, replicates=True, sample_count_frac: float = 0.3):

    if replicates:
        # adjust the sample counts to account replicates
        variant_sample_counts.update((x, y / 2) for x, y in variant_sample_counts.items())

    recurrent_variant_mutations = {}
    for variant, sample_count in variant_sample_counts.items():
        # get data for variant
        variant_df = df[df['VOC'] == variant]
        variant_isnvs = get_common_muts(variant_df, replicates)
        prominent_mutations = set(
            variant_isnvs[variant_isnvs['source_count'] > (sample_count_frac * sample_count)]['mutation'])
        recurrent_variant_mutations[variant] = prominent_mutations

    return recurrent_variant_mutations


def create_common_mutations_report(df, output_path, variant_mutations, file_title=None, sample_count_frac: float = 0.2,
                                   drop_gamma=True):
    mutations_report = pd.DataFrame(columns=['mutation', 'alpha', 'delta', 'gamma', 'omicron', 'non-voc', 'total'])

    # iterate over all the unique mutations in the data
    for mutation in set.union(*variant_mutations.values()):
        row = {'mutation': mutation, 'POS': df.loc[df['mutation'] == mutation, 'POS'].iloc[0],
               'REF': df.loc[df['mutation'] == mutation, 'REF'].iloc[0],
               'ALT': df.loc[df['mutation'] == mutation, 'ALT'].iloc[0],
               'GFF_FEATURE': df.loc[df['mutation'] == mutation, 'GFF_FEATURE'].iloc[0],
               'REF_AA': df.loc[df['mutation'] == mutation, 'REF_AA'].iloc[0],
               'ALT_AA': df.loc[df['mutation'] == mutation, 'ALT_AA'].iloc[0]}

        # get data about the mutation's aa
        row['POS_AA'], row['GFF_FEATURE'] = consensus_utils.nucleotide_to_aa_position(row['POS'],
                                                                                      gene_id=row['GFF_FEATURE'])

        # count the number of samples with the mutation for each variant
        for variant, mutations_list in variant_mutations.items():
            if mutation in mutations_list:
                count = df[(df['VOC'] == variant) & (df['mutation'] == mutation)]['sample'].nunique()
                row[variant.lower()] = count
            else:
                row[variant.lower()] = 0
        row['total'] = sum(row[variant.lower()] for variant in variant_mutations)

        # add synonymity column
        row['synonymity'] = masked_variants_utils.categorize_synonymity(row)

        # concat to report
        mutations_report = pd.concat([mutations_report, pd.DataFrame([row])], ignore_index=True)

    if drop_gamma:  # drop all mutations that are common only to gamma
        muts_to_drop = mutations_report[
            (mutations_report['alpha'] == 0) & (mutations_report['delta'] == 0) & (mutations_report['omicron'] == 0) & (
                        mutations_report['non-voc'] == 0) & (mutations_report['gamma'] > 0)]
        mutations_report.drop(muts_to_drop.index, inplace=True)

    # add CDS categorization
    mutations_report['CDS'] = mutations_report.apply(lambda row: masked_variants_utils.categorize_cds(row), axis=1)

    # add %GC in the 20 nucleotides surrounding the mutation
    mutations_report['GC_content_nearby'] = mutations_report.apply(
        lambda row: consensus_utils.get_GC_content_near_pos(row['POS'],
                                                            cons_seq_file='/Users/shir/Documents/University/Projects/AccuNGS_vs_iVar/MN908947_ref.fasta',
                                                            window_size=10), axis=1)

    # add distance from primers to the report
    mutations_report[['nt_from_V3_primers', 'V3_strand']] = mutations_report.apply(lambda row: pd.Series(
        consensus_utils.get_distance_from_closest_primer(row['POS'], 'V3',
                                                         '/Users/shir/Documents/University/Projects/acute_covid/primer_data')),
                                                                                   axis=1)
    mutations_report[['nt_from_V4.1_primers', 'V4.1_strand']] = mutations_report.apply(lambda row: pd.Series(
        consensus_utils.get_distance_from_closest_primer(row['POS'], 'V4.1',
                                                         '/Users/shir/Documents/University/Projects/acute_covid/primer_data')),
                                                                                       axis=1)

    # adjust the ORF1a & ORF1b positions to the ORF1ab gene
    mutations_report_fixed_orf1ab = masked_variants_utils.adjust_orf1ab_pos_aa_and_gff_feature(mutations_report)

    # add Bloom & Neher's fitness estimations to the report
    fitness_inference = pd.read_table(
        '/Users/shir/Documents/University/Projects/acute_covid/bloom_n_neher_fitness_inference/aamut_fitness_all.csv',
        sep=",")
    fitness_inference.rename(columns={'gene': 'GFF_FEATURE', 'mutant_aa': 'ALT_AA', 'aa_site': 'POS_AA'}, inplace=True)
    mutations_report_w_fitness = pd.merge(mutations_report_fixed_orf1ab,
                                          fitness_inference[['GFF_FEATURE', 'ALT_AA', 'POS_AA', 'delta_fitness']],
                                          how='left', on=['GFF_FEATURE', 'ALT_AA', 'POS_AA'])

    mutations_report_w_fitness.to_csv(
        f'{output_path}/{file_title}mutations_summary_{sample_count_frac}_sample_frac.tsv',
        sep="\t", index=False)
    return mutations_report_w_fitness


def main():
    # retrieve data & add necessary columns
    all_variants = pd.read_table('data/all_variants.tsv', sep="\t", index_col=0)
    utils.add_data_columns(all_variants)

    # merge with info files & get sample count of each variant before filtering
    all_variants_w_info = utils.merge_variants_w_info(all_variants,
                                                      sym_onset_path="data/sym_onset.csv",
                                                      ct_values_path="data/ct_values.csv")
    variant_sample_counts = all_variants_w_info.groupby('VOC')['sample'].nunique().to_dict()    # this count disregards replicates

    # get minor alleles
    minor_alleles = all_variants.loc[(all_variants['ALT_FREQ'] < 0.2)]

    # filter to eliminate as many sequencing errors as possible
    filtered_minor_alleles = filtering_functions.filter_mutations(minor_alleles, coverage_t=0, frequency_t=0.01,
                                                                  base_count_t=0, ignore_indels=False)

    # get suspicious recurrent mutations
    recurrent_mutations = find_recurrent_mutations_in_variants(filtered_minor_alleles, variant_sample_counts,
                                                               replicates=False, sample_count_frac = 0.3)
    # analyze only the recurrent mutations of variants that have at least 10 samples
    filtered_recurrent_mutations = {var: muts for var, muts in recurrent_mutations.items()
                                    if var in variant_sample_counts and variant_sample_counts[var] >= 10}
    # TODO: edit & implement the mutation report function
    mutations_report =
    pass


if __name__ == "__main__":
    main()
