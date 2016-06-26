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
    population_id_dict = {}
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
        population_id_dict.update({_[0]:j+1 for _ in new_keys})
        translation_dicts.append(translation_dict)
        inverse_translation_dict.update(dict(zip(translation_dict.values(), translation_dict.keys())))

    return anonymized_dict, inverse_translation_dict, population_id_dict


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


def write_geno_data(input_fn, geno_dict):
    df = pd.DataFrame(dict((k[0], hap1 + hap2) for k, (hap1, hap2) in geno_dict.items()))
    zfill_amt = int(np.log10(len(df)) - .9) + 1
    df.index = ['SNP_{{:0{}}}'.format(zfill_amt).format(_) for _ in df.index]
    df.to_csv(input_fn, sep=' ')

    # even_cols = [_ for _ in df.columns[::2]]
    # last_half_indivs = [_ for _ in df.index[len(df.index)/2:]]
    # masked_df = df.copy()
    # masked_df.ix[last_half_indivs, even_cols] = -1
    # masked_df.to_csv(input_fn, sep=' ')
    return


def write_population_data(input_fn, output_fn, population_dict, known_populations):
    df = pd.Series(population_dict)
    df.to_csv(output_fn, sep=' ')

    if known_populations:
        df.ix[len(df)/2:] = -1
        df.to_csv(input_fn, sep=' ')
    return


def make_admixed_individuals(indiv_dict, population_dict, n_indivs):
    n_segments = np.random.poisson(4, n_indivs)
    genome_length = len(indiv_dict.values()[0])
    # print n_segments[0]
    all_breakpoints = [[0] + sorted(np.random.choice(range(1, genome_length-1), n, replace=False)) + [genome_length] for n in n_segments]
    # print all_breakpoints[0]
    run_lengths = [[end - start for end, start in zip(breakpoints[1:], breakpoints[:-1])] for breakpoints in all_breakpoints]
    # print run_lengths[0]
    all_indivs_indices = [np.random.choice(range(len(indiv_dict)), n+1) for n in n_segments]
    # print all_indivs_indices[0]
    indiv_names = indiv_dict.keys()
    all_ancestors_keys = [[indiv_names[_][0] for _ in this_indiv_indices] for this_indiv_indices in all_indivs_indices]
    # print all_ancestors_keys[0]
    # print population_dict.items()[:10]
    all_populations = [[population_dict[_]
                        for _ in this_indiv_ancestor_keys] for this_indiv_ancestor_keys in all_ancestors_keys]
    # print all_populations[0]

    all_population_SNP_iterators = [itertools.chain.from_iterable([pop]*run_length for pop, run_length in zip(these_pops, these_run_lengths))
                                    for these_pops, these_run_lengths in zip(all_populations, run_lengths)]
    # print all_population_SNP_iterators[0]
    all_population_SNPs = [[_ for _ in SNP_iterator] for SNP_iterator in all_population_SNP_iterators]
    # print all_population_SNPs[0]

    starts = [_[:-1] for _ in all_breakpoints]
    # print starts[0]
    ends = [_[1:] for _ in all_breakpoints]
    # print ends[0]
    all_indivs_keys = [[indiv_names[_] for _ in this_indiv_indices] for this_indiv_indices in all_indivs_indices]
    # print all_indivs_keys[0]
    # print [indiv_dict[indiv] for indiv in all_indivs_keys[0]]
    z = [zip(_1, _2, _3) for _1, _2, _3 in zip(all_indivs_keys, starts, ends)]
    # print z[0]
    # print len(z)
    # print [indiv_dict[indiv] for indiv, start, end in z[0]]
    # print sum([len(indiv_dict[indiv][start:end]) for indiv, start, end in z[0]])

    all_SNPs_iterators = [itertools.chain.from_iterable([indiv_dict[indiv][start:end] for indiv, start, end in indivs_starts_ends])
                          for indivs_starts_ends in z]
    all_SNPs = [[_ for _ in SNP_iterator] for SNP_iterator in all_SNPs_iterators]
    all_SNPs = np.array(all_SNPs)
    random_noise = np.random.binomial(1, .01, size=all_SNPs.shape)
    all_SNPs += random_noise
    all_SNPs = all_SNPs % 3
    # print all_SNPs[0]
    # print len(all_SNPs[0])
    zfill_amt = int(np.log10(n_indivs))
    geno_dict = {'LOCAL_IMP_{{:0{}}}'.format(zfill_amt).format(i): (_1, _2) for _1, _2, i in zip(all_SNPs, all_population_SNPs, range(n_indivs))}
    return geno_dict


def write_admixed_indivs(input_fn, output_fn, admixed_dict):
    # print admixed_dict.items()[0]
    with open(input_fn, 'w') as f_in:
        with open(output_fn, 'w') as f_out:
            for k, v in sorted(admixed_dict.items()):
                f_in.write('{}:{}\n'.format(k, ''.join(str(_) for _ in v[0])))
                f_out.write('{}:{}\n'.format(k, ''.join(str(_) for _ in v[1])))
    return


def make_data(folder_name, n_snps, generation_sizes, additional_unrelated_indivs,
              n_populations, local, known_populations, n_admixed_indivs=0):
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

            # print len(parent_dict)
            # print len(unrelated_dict)
            population_dict = dict(parent_dict.items() + unrelated_dict.items())
            input_dicts.append(population_dict)

        complete_input_dict, inverse_translation_dict, population_ids = combine_populations(input_dicts)
        # print population_ids.items()[:10]

        geno_input_fn = join(folder_name, dataset, 'input', 'genotypes.txt')
        pop_input_fn = join(folder_name, dataset, 'input', 'populations.txt')
        output_fn = join(folder_name, dataset, 'output', 'population_tags.txt')

        write_geno_data(input_fn=geno_input_fn, geno_dict=complete_input_dict)
        write_population_data(input_fn=pop_input_fn, output_fn=output_fn,
                              population_dict=population_ids, known_populations=known_populations)

        if local:
            geno_dict = {k: v1 + v2 for k, (v1, v2) in complete_input_dict.items()}
            admixed_indivs = make_admixed_individuals(geno_dict, population_ids, n_admixed_indivs)
            admixed_input_fn = join(folder_name, dataset, 'input', 'admixed_genotypes.txt')
            admixed_output_fn = join(folder_name, dataset, 'output', 'admixed_populations.txt')
            write_admixed_indivs(input_fn=admixed_input_fn, output_fn=admixed_output_fn, admixed_dict=admixed_indivs)


if __name__ == "__main__":
    # Global Ancestry
    make_data(folder_name='easy', n_snps=1000, generation_sizes=(200,),
              additional_unrelated_indivs=0, n_populations=4, local=False, known_populations=True)

    make_data(folder_name='hard', n_snps=1000, generation_sizes=(200,),
              additional_unrelated_indivs=0, n_populations=8, local=False, known_populations=False)

    # Local Ancestry
    make_data(folder_name='medium', n_snps=10000, generation_sizes=(100,),
              additional_unrelated_indivs=0, n_populations=4, local=True, known_populations=True, n_admixed_indivs=100)

    make_data(folder_name='very_hard', n_snps=10000, generation_sizes=(100,),
              additional_unrelated_indivs=0, n_populations=8, local=True, known_populations=False, n_admixed_indivs=200)
