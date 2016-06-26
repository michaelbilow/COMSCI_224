import numpy as np
import itertools
from os.path import join, split, splitext, exists
from os import makedirs


def make_parents_and_unrelated(n_snps, n_genos, n_unrelated=0, population_index=0):
    no_break = False
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
                    no_break = True
        if not no_break: break
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
    zfill_amt = int(np.log10(len(ancestry_dict) + 1)) + 1
    translation_dict = {k[0]: 'IND_{}:0{}{}'.format('{', zfill_amt, '}').format(i) for k, i in
                        zip(ancestry_dict.keys(), range(1, len(ancestry_dict) + 1))}
    old_keys, vals = zip(*ancestry_dict.items())
    new_keys = [(translation_dict[indiv_id],
                 tuple(
                     tuple(translation_dict[ancestor_id] for ancestor_id in ancestry_level)
                     for ancestry_level in ancestry))
                for indiv_id, ancestry in old_keys]
    anonymized_dict = dict(zip(new_keys, vals))
    return anonymized_dict, translation_dict


def find_common_ancestors(ancestry1, ancestry2):
    all_ancestors1, all_ancestors2 = [_ for _ in itertools.chain.from_iterable(ancestry1)], \
                                     [_ for _ in itertools.chain.from_iterable(ancestry2)]
    common_ancestors = list(_ for _ in all_ancestors1 if _ in all_ancestors2)
    return common_ancestors


def is_child(child_id, parent_ancestry):
    common_ancestry_level = list(child_id in ancestry_level for ancestry_level in parent_ancestry)
    if any(common_ancestry_level):
        common_ancestry = common_ancestry_level.index(True) + 1
    else:
        common_ancestry = False
    return common_ancestry


def any_common_ancestry(indiv1, indiv2):
    id1, id2 = indiv1[0], indiv2[0]
    ancestry1, ancestry2 = indiv1[1], indiv2[1]
    return is_child(id1, ancestry2) or is_child(id2, ancestry1) or find_common_ancestors(ancestry1, ancestry2)


def is_sibling(ancestry1, ancestry2):
    parents1, parents2 = ancestry1[0], ancestry2[0]
    return parents1 and parents2 and all(_ in ancestry2[0] for _ in ancestry1[0])


def determine_relatedness(ancestry_dict):
    output_dict = {}
    for indiv1, indiv2 in itertools.combinations(ancestry_dict, 2):
        if indiv1 > indiv2:
            indiv1, indiv2 = indiv2, indiv1

        id1, id2 = indiv1[0], indiv2[0]
        ancestry1, ancestry2 = indiv1[1], indiv2[1]
        parent_status = is_child(id2, ancestry1)
        child_status = is_child(id1, ancestry2)
        if parent_status:
            relation = 'P{}'.format(parent_status)
        elif child_status:
            relation = 'C{}'.format(child_status)
        elif is_sibling(ancestry1, ancestry2):
            relation = 'S'
        else:
            common_ancestors = find_common_ancestors(ancestry1, ancestry2)
            if not common_ancestors:
                continue
            common_ancestry_level1 = [any([_ in ancestry_level for _ in common_ancestors]) for ancestry_level in
                                      ancestry1]
            common_ancestry_level2 = [any([_ in ancestry_level for _ in common_ancestors]) for ancestry_level in
                                      ancestry2]
            common_ancestry1 = common_ancestry_level1.index(True)
            common_ancestry2 = common_ancestry_level2.index(True)
            relation = 'M{}D{}'.format(common_ancestry1, common_ancestry2 - common_ancestry1)
        output_dict[(id1, id2)] = relation
    return output_dict


def combine_populations(population_dicts):
    total_n_individuals = sum(len(_) for _ in population_dicts)
    zfill_amt = int(np.log10(total_n_individuals - .9)) + 1
    translation_dicts = []
    anonymized_dict = {}
    inverse_translation_dict = {}
    print total_n_individuals
    new_keys = ['IND_{}:0{}{}'.format('{', zfill_amt, '}').format(i) for i in range(1, total_n_individuals+1)]
    shuffled_new_keys = np.random.permutation(new_keys)
    for j in range(len(population_dicts)):
        shuffled_new_keys_j = shuffled_new_keys[sum(len(_) for _ in population_dicts[:j]): sum(len(_) for _ in population_dicts[:j + 1])]
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


def make_data(folder_name, n_snps, generation_sizes, additional_unrelated_indivs, return_frontier=False,
              n_populations=1, random_fraction=1):
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

        if return_frontier:
            frontier_generation_size = max(len(_[1]) for _ in complete_input_dict.keys())
            complete_input_dict = {k: v for k, v in complete_input_dict.items()
                               if inverse_translation_dict[k[0]] < 0 or
                               len(k[1]) == frontier_generation_size}

        if random_fraction != 1:
            complete_input_dict_items = complete_input_dict.items()
            shuffled_complete_input_dict_items = np.random.permutation(range(len(complete_input_dict_items)))
            output_pop_size = int(random_fraction*len(shuffled_complete_input_dict_items))
            complete_input_dict = dict(complete_input_dict_items[_] for _ in shuffled_complete_input_dict_items[:output_pop_size])

        relatedness_dict = determine_relatedness(ancestry_dict=complete_input_dict)


        input_fn = join(folder_name, dataset, 'input', 'genotypes.txt')
        output_fn = join(folder_name, dataset, 'output', 'relationships.txt')


        write_input_data(fn=input_fn, ancestry_dict=complete_input_dict)
        write_output_data(fn=output_fn, relatedness_dict=relatedness_dict)


if __name__ == "__main__":
    make_data(folder_name='easy', n_snps=1000, generation_sizes=(100, 200),
              additional_unrelated_indivs=300, return_frontier=True)

    # make_data(folder_name='easy_all_bells', n_snps=1000, generation_sizes=(10, 20),
    #           additional_unrelated_indivs=30, return_frontier=True, n_populations=3, random_fraction=.7)

    make_data(folder_name='medium', n_snps=1000, generation_sizes=(100, 200, 300), additional_unrelated_indivs=400,
              return_frontier=False)

    make_data(folder_name='hard', n_snps=1000, generation_sizes=(100, 200, 300), additional_unrelated_indivs=400,
              return_frontier=False, n_populations=3)


    make_data(folder_name='very_hard', n_snps=1000, generation_sizes=(100, 200, 300), additional_unrelated_indivs=400,
              return_frontier=False, n_populations=6, random_fraction=.7)
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
