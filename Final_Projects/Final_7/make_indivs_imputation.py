import numpy as np
import itertools
from os.path import join, split, splitext, exists
from os import makedirs
import pandas as pd


def make_parents_and_unrelated(n_snps, n_genos, n_unrelated=0, population_index=0):
    identical_individuals = False
    n_indivs = n_genos + n_unrelated
    while True:
        freqs = np.random.uniform(low=(5.0 / n_snps), high=1 - (5.0 / n_snps), size=n_snps)
        haps = np.array([np.random.binomial(1, freq, 2 * n_indivs) for freq in freqs]).T
        for i in range(len(haps)):
            hap1 = haps[i]
            for j in range(i):
                hap2 = haps[j]
                # print i, j
                if (hap1 == hap2).all():
                    identical_individuals = True
        if not identical_individuals:
            break
    parent_haps = [(haps[2 * i], haps[2 * i + 1]) for i in range(n_genos)]
    parent_ancestors = [(i, ((),)) for i in range(n_genos)]
    ancestor_dict = dict(zip(parent_ancestors, parent_haps))

    unrelated_haps = [(haps[2 * i], haps[2 * i + 1]) for i in range(n_genos, n_indivs)]
    unrelated_ancestors = [(-i, ((),)) for i in range(1, n_unrelated + 1)]
    unrelated_dict = dict(zip(unrelated_ancestors, unrelated_haps))
    return ancestor_dict, unrelated_dict


def make_gamete(indiv, recombinations_lambda):
    starting_strand = np.random.binomial(1, .5)
    n_recombinations = np.random.poisson(recombinations_lambda)
    breakpoints = np.random.choice(range(1, len(indiv[0]) - 1), n_recombinations, replace=False).tolist()
    breakpoints = [0] + breakpoints + [len(indiv[0])]
    gamete_chunks = [indiv[strand][start:end] for strand, start, end in
                     zip(((starting_strand + i) % 2 for i in range(n_recombinations + 1)),
                         breakpoints[:-1],
                         breakpoints[1:])
                     ]
    gamete = np.array(list(itertools.chain.from_iterable(gamete_chunks)))
    return gamete


def make_child_genome(parent1, parent2, recombinations=0):
    gamete1 = make_gamete(parent1, recombinations)
    gamete2 = make_gamete(parent2, recombinations)
    return gamete1, gamete2


def make_child_ancestry(parent_pair):
    parent_ids = (parent_pair[0][0], parent_pair[1][0])
    child_ancestry = zip(parent_pair[0][1], parent_pair[1][1])
    # print list(_ for _ in child_ancestry)
    child_ancestry = [tuple(itertools.chain.from_iterable(_)) for _ in child_ancestry]
    # print child_ancestry
    child_ancestry = [parent_ids] + list(child_ancestry)

    child_ancestry = tuple(child_ancestry)
    # print child_ancestry
    return child_ancestry


def make_kids(ancestor_dict, n_kids, no_incest=True, monogamous=True, no_intergenerational=True):
    parents = ancestor_dict.keys()
    if no_intergenerational:
        frontier_generation_size = max(len(_[1]) for _ in parents)
        parents = [_ for _ in parents if len(_[1]) == frontier_generation_size]
    parent_pairs = list(itertools.combinations(parents, 2))
    if no_incest:
        parent_pairs = [(_1, _2) for _1, _2 in parent_pairs if not any_common_ancestry(_1, _2)]
        if not parent_pairs:
            raise ValueError
    if monogamous:
        parent_pairs = select_disjoint_pairs(parent_pairs)
        if not parent_pairs:
            raise ValueError
    breeding_pair_inds = np.random.choice(np.arange(len(parent_pairs)), n_kids, replace=True)
    breeding_pairs = [parent_pairs[_] for _ in breeding_pair_inds]
    kid_genomes = [make_child_genome(ancestor_dict[_1], ancestor_dict[_2], 0) for _1, _2 in breeding_pairs]
    kid_ids = [(len(ancestor_dict) + i, make_child_ancestry(pair)) for i, pair in
               zip(range(n_kids),
                   breeding_pairs)]
    for k, v in zip(kid_ids, kid_genomes):
        ancestor_dict[k] = v
    # print kid_ids
    # print kid_genomes
    return


def select_disjoint_pairs(pairs):
    all_elements = set()
    parent_pairs = []
    for pair in pairs:
        if any(_ in all_elements for _ in pair):
            continue
        else:
            parent_pairs.append(pair)
            all_elements.update(pair)
    return parent_pairs


def make_anonymized_genos_dict(ancestry_dict):
    zfill_amt = int(np.log10(len(ancestry_dict) + .1)) + 1
    translation_dict = {k[0]: 'IND_{{:0{}}}'.format(zfill_amt).format(i) for k, i in
                        zip(ancestry_dict.keys(), range(1, len(ancestry_dict) + 1))}
    old_keys, vals = zip(*ancestry_dict.items())
    new_keys = [(translation_dict[indiv_id],
                 tuple(
                     tuple(translation_dict[ancestor_id] for ancestor_id in ancestry_level)
                     for ancestry_level in ancestry))
                for indiv_id, ancestry in old_keys]
    anonymized_dict = dict(zip(new_keys, vals))
    return anonymized_dict, translation_dict


def combine_populations(population_dicts):
    total_n_individuals = sum(len(_) for _ in population_dicts)
    zfill_amt = int(np.log10(total_n_individuals + .1)) + 1
    translation_dicts = []
    anonymized_dict = {}
    inverse_translation_dict = {}
    print total_n_individuals
    new_keys = ['IND_{{:0{}}}'.format(zfill_amt).format(i) for i in range(1, total_n_individuals + 1)]
    shuffled_new_keys = np.random.permutation(new_keys)
    for j in range(len(population_dicts)):
        shuffled_new_keys_j = shuffled_new_keys[
                              sum(len(_) for _ in population_dicts[:j]): sum(len(_) for _ in population_dicts[:j + 1])]
        translation_dict = dict(zip((_[0] for _ in population_dicts[j].keys()), shuffled_new_keys_j))
        print len(translation_dict)
        old_keys, vals = zip(*population_dicts[j].items())
        new_keys = [(translation_dict[indiv_id],
                     tuple(
                         tuple(translation_dict[ancestor_id] for ancestor_id in ancestry_level)
                         for ancestry_level in ancestry))
                    for indiv_id, ancestry in old_keys]
        anonymized_dict_update = dict(zip(new_keys, vals))
        anonymized_dict.update(anonymized_dict_update)
        translation_dicts.append(translation_dict)
        inverse_translation_dict.update(dict(zip(translation_dict.values(), translation_dict.keys())))

    return anonymized_dict, inverse_translation_dict


def make_phenos(genotypes_dict, N_snp_effects, N_snp_snp_effects):
    sorted_indivs = sorted(genotypes_dict.keys())
    sorted_genos = [genotypes_dict[_] for _ in sorted_indivs]

    genotypes = np.array([hap1 + hap2 for hap1, hap2 in sorted_genos])
    phenotypes = np.random.normal(size=genotypes.shape[0])
    normalized_genotypes = (genotypes - np.mean(genotypes, axis=0)) / np.std(genotypes, axis=0)
    snp_effect_sizes = np.linspace(.1, 3, N_snp_effects)
    snp_snp_effect_sizes = np.linspace(.1, 5, N_snp_snp_effects)

    SNPs = range(genotypes.shape[1])
    SNP_SNPs = [_ for _ in itertools.combinations(SNPs, 2)]

    snp_effect_ids = np.random.permutation(SNPs)[:N_snp_effects]
    snp_snp_effect_ids = np.random.permutation(SNP_SNPs)[:N_snp_snp_effects]
    normalized_snp_values = normalized_genotypes[:, snp_effect_ids]
    snp_snp_values = np.array(
        [normalized_genotypes[:, _1] * normalized_genotypes[:, _2] for _1, _2 in snp_snp_effect_ids]).T
    normalized_snp_snp_values = (snp_snp_values - np.mean(snp_snp_values, axis=0)) / np.std(snp_snp_values, axis=0)

    snp_effects = np.dot(normalized_snp_values, snp_effect_sizes)
    snp_snp_effects = np.dot(normalized_snp_snp_values, snp_snp_effect_sizes)
    print snp_effects
    print snp_snp_effects

    phenotypes += snp_effects + snp_snp_effects
    pheno_dict = {k[0]: v for k, v in zip(sorted_indivs, phenotypes)}

    print snp_effect_ids
    print snp_effect_sizes
    print snp_snp_effect_ids
    print snp_snp_effect_sizes
    associations_dict = dict(itertools.chain(zip(snp_effect_ids,
                                                 snp_effect_sizes),
                                             zip(((_1, _2) for _1, _2 in snp_snp_effect_ids),
                                                 snp_snp_effect_sizes)))
    return pheno_dict, associations_dict


def write_input_data(fn, ancestry_dict):
    # print ancestry_dict
    with open(fn, 'w') as f:
        for (iid, ancestry), (haplotype1, haplotype2) in sorted(ancestry_dict.items()):
            f.write('{}:{}\n'.format(iid, ''.join((haplotype1 + haplotype2).astype(str))))


def write_output_data(fn, relatedness_dict):
    with open(fn, 'w') as f:
        sorted_relations = sorted(relatedness_dict.keys())
        for relation in sorted_relations:
            f.write('{}-{}:{}\n'.format(relation[0], relation[1], relatedness_dict[relation]))


def write_geno_data(input_fn, output_fn, geno_dict):
    df = pd.DataFrame(dict((k[0], hap1 + hap2) for k, (hap1, hap2) in geno_dict.items()))
    zfill_amt = int(np.log10(len(df)) - .9) + 1
    df.index = ['SNP_{{:0{}}}'.format(zfill_amt).format(_) for _ in df.index]
    df.to_csv(output_fn, sep=' ')

    even_cols = [_ for _ in df.columns[::2]]
    last_half_indivs = [_ for _ in df.index[len(df.index)/2:]]
    masked_df = df.copy()
    masked_df.ix[last_half_indivs, even_cols] = -1
    masked_df.to_csv(input_fn, sep=' ')
    return




def make_data(folder_name, n_snps, generation_sizes, additional_unrelated_indivs,
              n_populations):
    all_folders = [join(folder_name, dataset, input_or_output)
                   for dataset in ('example', 'test') for input_or_output in ('input', 'output')]
    for _ in all_folders:
        if not exists(_):
            makedirs(_)
    for dataset in ('example', 'test'):
        input_dicts = []
        for pop in range(n_populations):
            parent_dict, unrelated_dict = make_parents_and_unrelated(n_snps=n_snps,
                                                                     n_genos=generation_sizes[0],
                                                                     n_unrelated=additional_unrelated_indivs)

            for generation_size in generation_sizes[1:]:
                make_kids(ancestor_dict=parent_dict, n_kids=generation_size)

            print len(parent_dict)
            print len(unrelated_dict)
            population_dict = dict(parent_dict.items() + unrelated_dict.items())
            input_dicts.append(population_dict)

        complete_input_dict, inverse_translation_dict = combine_populations(input_dicts)

        input_fn = join(folder_name, dataset, 'input', 'raw_genotypes.txt')
        output_fn = join(folder_name, dataset, 'output', 'imputed_genotypes.txt')

        write_geno_data(input_fn=input_fn, output_fn=output_fn, geno_dict=complete_input_dict)


if __name__ == "__main__":
    make_data(folder_name='easy', n_snps=100, generation_sizes=(500,),
              additional_unrelated_indivs=0, n_populations=1)

    make_data(folder_name='medium', n_snps=10000, generation_sizes=(500,),
              additional_unrelated_indivs=0, n_populations=1)

    make_data(folder_name='hard', n_snps=10000, generation_sizes=(250,),
              additional_unrelated_indivs=0, n_populations=4)

    # make_data(folder_name='easy_all_bells', n_snps=1000, generation_sizes=(10, 20),
    #           additional_unrelated_indivs=30, return_frontier=True, n_populations=3, random_fraction=.7)

    # make_data(folder_name='medium', n_snps=1000, generation_sizes=(100, 200, 300), additional_unrelated_indivs=400,
    #           return_frontier=False)
    #
    # make_data(folder_name='hard', n_snps=1000, generation_sizes=(100, 200, 300), additional_unrelated_indivs=400,
    #           return_frontier=False, n_populations=3)


    # make_data(folder_name='very_hard', n_snps=1000, generation_sizes=(100, 200, 300), additional_unrelated_indivs=400,
    #           return_frontier=False, n_populations=6, random_fraction=.7)
    # d, unrelated_d = make_parents_and_unrelated(n_snps=500, n_genos=100)
    # # for k in d:
    # #     print k
    # #     print d[k]
    # make_kids(d, 100)
    #
    # # for k in d:
    # #     print k
    # #     print d[k]
    #
    # make_kids(d, 100)
    # # for k, v in d.items():
    # #     print k
    # #     print v
    #
    # new_d = make_anonymized_genos_dict(d)
    # # for k in new_d:
    # #     print k
    #
    # relatedness_d = determine_relatedness(new_d)
    # for k, v in relatedness_d.items():
    #     print k, v
    #
    # write_input_data('test.txt', new_d)
    # write_output_data('test_output.txt', relatedness_d)
