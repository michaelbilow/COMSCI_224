import pandas as pd
import numpy as np
from scipy.stats import uniform, norm
from scipy.sparse import block_diag


def make_individuals():
    raw_indivs = norm.rvs(size=(N_INDIVIDUALS, N_SNPS))
    if CORR_MAT_SQRT:
        raw_indivs = np.dot(CORR_MAT_SQRT, raw_indivs)
    else:
        pass
    one_or_more_minor_alleles = raw_indivs > SNP_THRESHOLDS[:, 0]
    # print one_or_more_minor_alleles
    two_minor_alleles = raw_indivs > SNP_THRESHOLDS[:, 1]
    genetic_state = one_or_more_minor_alleles.astype(int) + two_minor_alleles.astype(int)
    # print genetic_state
    # print genetic_state.mean(axis=0)
    # print genetic_state.mean(axis=0)/2 - MAFS
    # if abs((genetic_state.mean(axis=0)/2 - MAFS).mean()) > 3*(genetic_state.mean(axis=0)/2 - MAFS).std():
    #     raise ValueError

    zeroed_genetic_state = genetic_state/2 - MAFS
    # print zeroed_genetic_state
    # print EFFECT_SIZES

    liabilities = np.dot(zeroed_genetic_state, EFFECT_SIZES) + norm.rvs(size=indivs_per_round)

    all_cases = np.concatenate(cases)
    # print len(all_cases), n_cases
    all_controls = np.concatenate(controls)

    choose_cases = np.random.choice(len(all_cases), size=n_cases, replace=False)
    choose_controls = np.random.choice(len(all_controls), size=n_controls, replace=False)

    output_cases = all_cases[choose_cases, :]
    output_controls = all_controls[choose_controls]

    # print output_cases
    # print output_controls

    statuses = np.array(['Case']*len(output_cases) + ['Control'] * len(output_controls))

    output_df = pd.DataFrame(np.concatenate((output_cases, output_controls)))
    output_df['Status'] = statuses
    output_df = output_df.reindex(np.random.permutation(output_df.index))
    output_df.index = ['IND{:04}'.format(_) for _ in range(len(output_df))]
    output_df.columns = ['SNP{:05}'.format(_) for _ in range(len(output_df.columns)-1)] + ['Status']
    output_df.to_csv('SNP_Status.txt', sep=' ')
    import zipfile
    with zipfile.ZipFile('SNP_Status.zip', 'w') as myzip:
        myzip.write('SNP_status.txt')
    print output_df.head()
    return output_cases, output_controls



def make_corr_mat(corr_size):
    corr_mat = norm.rvs(size=(corr_size, corr_size))
    row_means = corr_mat.mean(axis=1)
    row_stds = corr_mat.std(axis=1)
    corr_mat = (corr_mat - row_means[:, np.newaxis]) / row_stds[:, np.newaxis]

    corr_mat = np.dot(corr_mat.T, corr_mat) / corr_size
    print corr_mat == corr_mat.T
    print corr_mat.shape
    print corr_mat.max()
    return

def make_weak_corr_mat(corr_size):
    raw_corr_mats = [np.ones((corr_size, corr_size)) * float(i)/10 for i in range(10)]
    corr_mats = [np.maximum(_, np.eye(corr_size)) for _ in raw_corr_mats]
    corr_mat_sqrts = [np.linalg.cholesky(_) for _ in corr_mats]
    output_corr_mats = np.random.choice(len(corr_mat_sqrts), N_SNPS/corr_size, replace=True)
    output_corr_mat_sqrts = [corr_mat_sqrts[_] for _ in output_corr_mats]
    output_corr_mat_sqrt = block_diag(output_corr_mat_sqrts)
    return output_corr_mat_sqrt


def make_populations(n_populations=6):
    import random

    def constrained_sum_sample_pos(n, total):
        """Return a randomly chosen list of n positive integers summing to total.
        Each such list is equally likely to occur."""

        dividers = sorted(random.sample(xrange(1, total), n - 1))
        return [a - b for a, b in zip(dividers + [total], [0] + dividers)]
    while True:
        corr_sizes = constrained_sum_sample_pos(n_populations, total=N_INDIVIDUALS)
        if all(_ > 100 for _ in corr_sizes):
            break
    raw_corr_mats = [np.ones((corr_size, corr_size)) * .6 for corr_size in corr_sizes]
    corr_mats = [np.maximum(_, np.eye(corr_size)) for _ in raw_corr_mats]
    corr_mat_sqrts = [np.linalg.cholesky(_) for _ in corr_mats]
    output_corr_mats = np.random.choice(len(corr_mat_sqrts), N_SNPS/corr_size, replace=True)
    output_corr_mat_sqrts = [corr_mat_sqrts[_] for _ in output_corr_mats]
    output_corr_mat_sqrt = block_diag(output_corr_mat_sqrts)
    return output_corr_mat_sqrt



def print_params():
    gl_dict = globals().copy()
    with open('params.txt', 'w') as f:
        for k in gl_dict:
            if k == k.upper():
                f.write('{}:{}\n'.format(k, gl_dict[k]))
    return

if __name__ == "__main__":
    N_INDIVIDUALS = 2000
    DISEASE_FREQUENCY = .1
    N_SNPS = 10000
    MAKE_CORRS = False

    CORRELATION_SIZE = 100

    if MAKE_CORRS:
        CORR_MAT_SQRT = make_weak_corr_mat(CORRELATION_SIZE)
    else:
        CORR_MAT_SQRT = False



    # assert N_SNPS % CORRELATION_SIZE == 0
    # FULL_CORR_MAT = np.kron(CORR_MAT, np.eye(N_SNPS/CORRELATION_SIZE))

    MAFS = uniform.rvs(0.05, .5, N_SNPS)
    print MAFS[:10]
    SNP_THRESHOLDS = np.array([np.array([norm.ppf((1.0 - maf)**2), norm.ppf(1.0 - maf**2)]) for maf in MAFS])
    print SNP_THRESHOLDS[:10]

    EFFECT_SIZES = np.random.permutation(np.concatenate((norm.rvs(0, .6, N_SNPS/1000), np.zeros(999*N_SNPS/1000))))
    # print_params()
    with open('true_effect_sizes.txt', 'w') as f:
        for _ in EFFECT_SIZES:
            f.write('{}\n'.format(_))
    make_individuals()