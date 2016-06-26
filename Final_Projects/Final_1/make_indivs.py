import numpy as np
import itertools
from os.path import join, split, splitext, exists
from os import makedirs


def make_haps(n_snps, n_haps, fraction=.5):
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


    # print haps
    hap_pairs = itertools.combinations_with_replacement(haps, 2)
    hap_pairs = [_ for _ in hap_pairs][::int(1.0/fraction)]
    # for _ in hap_pairs: print _
    indivs = [''.join((x + y).astype(str)) for x, y in hap_pairs]
    np.random.shuffle(indivs)
    return haps, indivs


def write_haplotypes(fn, haplotypes):
    with open(fn, 'w') as f:
        for hap in haplotypes:
            hap_str = ''.join(hap.astype(str))
            f.write('{}\n'.format(hap_str))
    return


def write_reads(fn, reads):
    with open(fn, 'w') as f:
        for read in reads:
            f.write('{}\n'.format(read))
    return


def make_data(folder_name,  n_snps, n_haps, fraction=.5, other_name=''):
    all_folders = [join(folder_name, dataset, input_or_output)
                   for dataset in ('example', 'test') for input_or_output in ('input', 'output')]
    for _ in all_folders:
        if not exists(_):
            makedirs(_)
    print 'Making {}{}'.format(folder_name, other_name if not other_name else ' ' + other_name)
    for dataset in ('example', 'test'):
        haplotypes, reads = make_haps(n_snps=n_snps, n_haps=n_haps, fraction=fraction)

        hap_fn = join(folder_name, dataset, 'output', 'haplotypes{}.txt'.format(other_name if not other_name else '_' + other_name))
        genotype_fn = join(folder_name, dataset, 'input', 'genotypes{}.txt'.format(other_name if not other_name else '_' + other_name))
        write_haplotypes(fn=hap_fn, haplotypes=haplotypes)
        write_reads(fn=genotype_fn, reads=reads)
    return


if __name__ == "__main__":
    make_data(folder_name='easy', n_snps=15, n_haps=6, fraction=.75)
    make_data(folder_name='medium', other_name='10_SNPs', n_snps=10, n_haps=6, fraction=.5)
    make_data(folder_name='medium', other_name='30_SNPs', n_snps=30, n_haps=20, fraction=.5)
    make_data(folder_name='medium', other_name='100_SNPs', n_snps=100, n_haps=25, fraction=.5)
    make_data(folder_name='medium', other_name='300_SNPs', n_snps=300, n_haps=30, fraction=.5)
    make_data(folder_name='hard', n_snps=1000, n_haps=40, fraction=.33)
    make_data(folder_name='very_hard', n_snps=10000, n_haps=80, fraction=.25)
