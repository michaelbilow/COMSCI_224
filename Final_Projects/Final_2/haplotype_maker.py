import pandas as pd
import numpy as np
from os.path import exists, join, split, splitext
from os import makedirs

SNP_FREQ_LEVELS = 5
MIN_READ_LENGTH = 3
READ_LENGTH = 50

def make_reads(genome_length, n_reads, avg_snps_per_read, n_chromosomes, error_rate):
    assert genome_length % SNP_FREQ_LEVELS == 0
    if n_chromosomes == 2:
        haplotype1 = np.random.binomial(1, .5, genome_length)
        haplotype2 = 1 - haplotype1
        haplotypes = np.array([haplotype1, haplotype2]).T
    else:
        haplotype_snps = [np.random.binomial(1, (1.0 + i) / (SNP_FREQ_LEVELS + 2.0),
                                             size=(genome_length / SNP_FREQ_LEVELS, n_chromosomes - 1))
                          for i in range(SNP_FREQ_LEVELS)]
        haplotypes = np.concatenate(haplotype_snps)
        np.random.shuffle(haplotypes)
        last_haplotype = 1 - haplotypes.max(axis=1)
        haplotypes = np.concatenate((haplotypes, last_haplotype.reshape(len(last_haplotype), 1)), axis=1)
        haplotypes = haplotypes.T
        np.random.shuffle(haplotypes)
        haplotypes = haplotypes.T
        np.random.shuffle(haplotypes)


    snp_dists = np.random.poisson(lam=float(READ_LENGTH) / avg_snps_per_read, size=n_reads - 1) + 1
    snp_locs = snp_dists.cumsum()
    snp_locs = np.array([0] + snp_locs.tolist())
    # print snp_locs
    read_lengths = [((snp_locs >= snp_loc) & (snp_locs < snp_loc + READ_LENGTH)).sum() for snp_loc in snp_locs]
    # print read_lengths
    # print np.mean(read_lengths)

    read_starts = np.random.randint(0, genome_length, size=n_reads)
    read_starts = sorted(read_starts)
    # read_lengths = np.random.poisson(lam=avg_read_length - MIN_READ_LENGTH, size=n_reads) + MIN_READ_LENGTH
    read_strands = np.random.randint(0, n_chromosomes, size=n_reads)
    reads = [haplotypes[read_starts[i]: read_starts[i] + read_lengths[read_starts[i]], read_strands[i]]
             for i in range(n_reads)]
    if error_rate != 0:
        reads = [(_ + np.random.binomial(1, error_rate, size=len(_))) % 2 for _ in reads]
    pretty_reads = ['{}{}{}'.format('-'*read_starts[i], ''.join(str(_) for _ in reads[i]),
                                    '-' * (genome_length - read_starts[i] - read_lengths[read_starts[i]]))
                    for i in range(n_reads)]

    # print reads[:10]
    # print read_starts[:10]
    # print all_pretty_reads[:20]
    return haplotypes, pretty_reads


def write_haplotypes(fn, haplotypes):
    hap_strs = (''.join(_) for _ in haplotypes.T.astype(str))
    with open(fn, 'w') as f:
        for _ in hap_strs:
            f.write('{}\n'.format(_))
    return


def write_reads(fn, reads):
    with open(fn, 'w') as f:
        for _ in reads:
            f.write('{}\n'.format(_))
    return


def create_data(folder_name, genome_length, n_reads, avg_snps_per_read, n_chromosomes=2, error_rate=0.0, other_name=''):
    all_folders = [join(folder_name, dataset, input_or_output)
                   for dataset in ('example', 'test') for input_or_output in ('input', 'output')]
    for _ in all_folders:
        if not exists(_):
            makedirs(_)
    print 'Making {}{}'.format(folder_name, other_name if not other_name else ' ' + other_name)
    for dataset in ('example', 'test'):
        haplotypes, reads = make_reads(genome_length=genome_length,
                                       n_reads=n_reads,
                                       avg_snps_per_read=avg_snps_per_read,
                                       n_chromosomes=n_chromosomes,
                                       error_rate=error_rate)

        input_fn = join(folder_name, dataset, 'input', 'reads{}.txt'.format(other_name if not other_name else '_' + other_name))
        output_fn = join(folder_name, dataset, 'output', 'haplotypes{}.txt'.format(other_name if not other_name else '_' + other_name))

        write_reads(fn=input_fn, reads=reads)
        write_haplotypes(fn=output_fn, haplotypes=haplotypes)
    return


if __name__ == "__main__":
    create_data(folder_name='easy', genome_length=100, n_reads=600,
                avg_snps_per_read=5, error_rate=0.0)
    create_data(folder_name='medium', other_name='low_error', genome_length=1000,
                n_reads=6000, avg_snps_per_read=5, error_rate=0.02)
    create_data(folder_name='medium', other_name='high_error', genome_length=1000,
                n_reads=6000, avg_snps_per_read=5, error_rate=0.1)
    create_data(folder_name='hard', other_name='4_chromosome', genome_length=1000,
                n_reads=10000, avg_snps_per_read=5, n_chromosomes=4, error_rate=0.02)
    create_data(folder_name='hard', other_name='8_chromosome', genome_length=1000,
                n_reads=18000, avg_snps_per_read=5, n_chromosomes=4, error_rate=0.02)
