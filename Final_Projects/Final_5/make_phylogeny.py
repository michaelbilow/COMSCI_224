import numpy as np
import itertools
import string
from os.path import join, split, splitext, exists
from os import makedirs


def make_haps(n_snps, n_haps, hap_names):
    no_break = False
    while True:
        freqs = np.random.uniform(low=(2.0/n_snps), high=1-(2.0/n_snps), size=n_snps)
        haps = np.array([np.random.binomial(1, freq, n_haps) for freq in freqs]).T
        for i in range(len(haps)):
            hap1 = haps[i]
            for j in range(i):
                hap2 = haps[j]
                # print i, j
                if (hap1 == hap2).all():
                    no_break = True
        if not no_break: break

    output_dict = {hap_names[i]: haps[i] for i in range(len(haps))}
    ancestry = {k: '' for k in output_dict}
    return output_dict, ancestry


def make_new_generation(pvs_haps, N_new_haps, hap_names, recombinations=0, mutation_rate=.05):
    N_no_recombination_haps = N_new_haps - recombinations
    base_on_hap_inds = np.random.choice(pvs_haps.keys(), N_no_recombination_haps)
    base_on_haps = [pvs_haps[ind] for ind in base_on_hap_inds]
    ancestral_haplotypes = base_on_hap_inds.tolist()

    recombination_hap_inds = np.random.choice(pvs_haps.keys(), size=(recombinations, 2), replace=False)
    recombined_haps = []
    hap_length = len(pvs_haps.values()[0])
    for i in range(recombinations):
        ancestral_hap1, ancestral_hap2 = recombination_hap_inds[i, :]
        hap1 = pvs_haps[ancestral_hap1]
        hap2 = pvs_haps[ancestral_hap1]
        breakpoint = np.random.randint(hap_length)
        recombined_hap = np.append(hap1[:breakpoint], hap2[breakpoint:])
        recombined_haps.append(recombined_hap)
        ancestral_haplotypes.append((ancestral_hap1, ancestral_hap2, breakpoint))
    if recombined_haps:
        base_on_haps = np.vstack((base_on_haps, recombined_haps))
    else:
        base_on_haps = np.array(base_on_haps)
    # ancestral_haplotypes = base_on_haps
    # print base_on_haps.shape

    n_bases_to_change = np.random.poisson(lam=mutation_rate*hap_length, size=N_new_haps) + 1
    bases_to_change = [np.random.choice(np.arange(hap_length), size=n, replace=False)
                       for n in n_bases_to_change]
    changes = [[1 if i in base_list else 0 for i in range(hap_length)]
               for base_list in bases_to_change]
    changes = np.array(changes)
    output_haps = (base_on_haps + changes) % 2

    output_dict = {hap_names[i]: output_haps[i] for i in range(len(output_haps))}
    new_ancestry = dict(zip(hap_names[:len(output_haps)],
                            ancestral_haplotypes))
    return output_dict, new_ancestry


def write_ancestry(ancestry_fn, haplotype_fn, ancestry, haplotype_dict):
    masking_dict = dict(zip(ancestry.keys(),
                            ('HAP_{:03}'.format(i)
                             for i in np.random.permutation(range(1, len(ancestry) + 1)))))
    # print masking_dict

    ancestry = {masking_dict[k]: '' if v == '' else masking_dict[v] if v in masking_dict else (masking_dict[v[0]], masking_dict[v[1]], v[2])
                for k, v in ancestry.items()}
    haplotype_dict = {masking_dict[k]: v for k,v in haplotype_dict.items()}

    with open(haplotype_fn, 'w') as f:
        for k in sorted(haplotype_dict.keys())[:2*len(haplotype_dict)/3]:
            hap_str = ''.join(haplotype_dict[k].astype(str))
            f.write('{}:{}\n'.format(k, hap_str))

    with open(ancestry_fn, 'w') as f:
        for k in sorted(haplotype_dict.keys()):
            f.write('{}:{}\n'.format(k, ancestry[k] if isinstance(ancestry[k], basestring)
            else ','.join(str(_) for _ in ancestry[k])))
    return


def make_data(folder_name, n_snps, n_haps, n_generations, recombinations_per_generation=0, mutation_rate=.05, other_name=''):
    all_folders = [join(folder_name, dataset, input_or_output)
                   for dataset in ('example', 'test') for input_or_output in ('input', 'output')]
    for _ in all_folders:
        if not exists(_):
            makedirs(_)

    hap_names = [''.join(_) for _ in itertools.product(string.uppercase, string.uppercase, string.uppercase)]
    print 'Making {}'.format(folder_name)
    for dataset in ('example', 'test'):
        pvs_haps, ancestry = make_haps(n_snps=n_snps, n_haps=n_haps, hap_names=hap_names)
        hap_names = hap_names[len(pvs_haps):]
        for i in range(2, n_generations + 1):
            print 'Generation {}'.format(i)
            new_haps, new_ancestry = make_new_generation(pvs_haps=pvs_haps,
                                                         N_new_haps=np.random.randint((2**i)*n_haps, 2**(i+1)*n_haps),
                                                         recombinations=recombinations_per_generation,
                                                         hap_names=hap_names,
                                                         mutation_rate=mutation_rate)
            hap_names = hap_names[len(new_haps):]
            ancestry = dict(ancestry.items() + new_ancestry.items())
            # print ancestry
            pvs_haps = dict(pvs_haps.items() + new_haps.items())
            # print pvs_haps

        hap_fn = join(folder_name, dataset, 'input', 'haplotypes{}.txt'.format(other_name if not other_name else '_' + other_name))
        ancestry_fn = join(folder_name, dataset, 'output', 'ancestry{}.txt'.format(other_name if not other_name else '_' + other_name))
        write_ancestry(ancestry_fn, hap_fn, ancestry, haplotype_dict=pvs_haps)


if __name__ == "__main__":
    make_data('easy', n_snps=100, n_haps=5, n_generations=3)
    # make_data('medium', n_snps=1000, n_haps=10, n_generations=4, recombinations_per_generation=2)
    # make_data('very_hard', n_snps=10000, n_haps=15, n_generations=5, recombinations_per_generation=3)
    # make_data('three_generations', n_snps=100, n_haps=3, n_generations=3)
    # make_data('two_gens_recombination', n_snps=100, n_haps=3, n_generations=2, recombinations_per_generation=1)
    # make_data('three_gens_recombination', n_snps=100, n_haps=4, n_generations=3, recombinations_per_generation=2)
    # make_data('hard', n_snps=1000, n_haps=20, n_generations=4, recombinations_per_generation=5)

