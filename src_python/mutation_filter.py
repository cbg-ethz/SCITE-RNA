"""
This code calculates the posterior probability of different mutation types and genotypes given the data
and can filter SNVs based on this posterior. Additionally, it allows to fit the parameters of the model
using the observed read counts and inferred genotypes.
"""

import numpy as np
from numba import njit
import math
from scipy.special import loggamma, logsumexp
from scipy.optimize import minimize
from scipy.stats import gamma, beta

from src_python.utils import load_config_and_set_random_seed

config = load_config_and_set_random_seed()


class MutationFilter:
    def __init__(self, error_rate=0.05, overdispersion=10, genotype_freq=None, mut_freq=0.5,
                 dropout_alpha=2, dropout_beta=8, dropout_direction_prob=0.5, overdispersion_h=6):
        """
       Filter mutations in single-cell RNA sequencing data and calculate the posterior probability of different
       mutation types given the observed alternative and reference read counts.

       Attributes:
           error_rate (float): The error rate for sequencing.
           overdispersion (float): The overdispersion parameter for homozygous genotypes.
           genotype_freq (dict): The frequency of each genotype.
           mut_freq (float): The mutation frequency.
           dropout_alpha (float): Alpha parameter for dropout probability.
           dropout_beta (float): Beta parameter for dropout probability.
           dropout_direction_prob (float): # Frequency, that the reference allele is dropped out instead of the alternative in the heterozygous case.
           overdispersion_h (float): The overdispersion parameter for heterozygous genotypes.
       """

        self.mut_type_prior = {s: np.nan for s in ['R', 'H', 'A', 'RH', 'HR', 'AH', 'HA']}
        if genotype_freq is None:
            genotype_freq = {'R': 1 / 3, 'H': 1 / 3, 'A': 1 / 3}
        self.genotype_freq = genotype_freq
        self.alpha_R = error_rate * overdispersion
        self.beta_R = overdispersion - self.alpha_R
        self.alpha_A = (1 - error_rate) * overdispersion
        self.beta_A = overdispersion - self.alpha_A
        self.alpha_H = 0.5 * overdispersion_h  # centered at 0.5, as the cells are assumed to be independent
        self.beta_H = overdispersion_h - self.alpha_H

        self.set_mut_type_prior(genotype_freq, mut_freq)
        self.dropout_prob = dropout_alpha / (dropout_alpha + dropout_beta)
        self.dropout_direction_prob = dropout_direction_prob
        self.overdispersion_H = overdispersion_h
        self.overdispersion = overdispersion
        self.error_rate = error_rate

    def update_alpha_beta(self, error_r, overdispersion_hom):
        self.alpha_R = error_r * overdispersion_hom
        self.beta_R = overdispersion_hom - self.alpha_R

        self.alpha_A = (1 - error_r) * overdispersion_hom
        self.beta_A = overdispersion_hom - self.alpha_A

    def set_mut_type_prior(self, genotype_freq, mut_freq):
        """
        Calculates and stores the log-prior for each possible mutation type of a locus (including non-mutated)

        [Arguments]
            genotype_freq: priors of the root (wildtype) having genotype R, H or A
            mut_freq: a priori proportion of loci that are mutated
        """
        # three non-mutated cases
        self.mut_type_prior['R'] = genotype_freq['R'] * (1 - mut_freq)
        self.mut_type_prior['H'] = genotype_freq['H'] * (1 - mut_freq)
        self.mut_type_prior['A'] = genotype_freq['A'] * (1 - mut_freq)
        # four mutated cases
        self.mut_type_prior['RH'] = genotype_freq['R'] * mut_freq
        self.mut_type_prior['HA'] = genotype_freq['H'] * mut_freq / 2  # can either become alternative or reference
        self.mut_type_prior['HR'] = self.mut_type_prior['HA']
        self.mut_type_prior['AH'] = genotype_freq['A'] * mut_freq

        # convert to log scale
        for s in self.mut_type_prior:
            self.mut_type_prior[s] = np.log(self.mut_type_prior[s])


    def single_read_llh_with_dropout(self, n_alt, n_total, genotype):
        """
        [Arguments]
            n_ref: number of ref reads
            n_total: total number of reads (ref + alt)
            genotype: the genotype of interest

        [Returns]
            the log-likelihood of observing n_ref, n_alt, given genotype
        """
        if genotype == 'R':
            result = log_betabinom_pmf(n_alt, n_total, self.alpha_R, self.beta_R)
        elif genotype == 'A':
            result = log_betabinom_pmf(n_alt, n_total, self.alpha_A, self.beta_A)
        elif genotype == 'H':
            result_no_dropout = np.log(1 - self.dropout_prob) + \
                                log_betabinom_pmf(n_alt, n_total, self.alpha_H, self.beta_H)
            dropout_R = np.log(self.dropout_prob) + np.log(1 - self.dropout_direction_prob) + \
                        log_betabinom_pmf(n_alt, n_total, self.alpha_R, self.beta_R)

            dropout_A = np.log(self.dropout_prob) + np.log(self.dropout_direction_prob) + \
                        log_betabinom_pmf(n_alt, n_total, self.alpha_A, self.beta_A)
            result = logsumexp([result_no_dropout, dropout_R, dropout_A])
        else:
            raise ValueError('[MutationFilter.single_read_llh] Invalid genotype.')

        return result


    def single_read_llh_with_individual_dropout(self, n_alt, n_total, genotype, dropout_prob,
                                                alpha_h, beta_h):
        """
        [Arguments]
            n_ref: number of ref reads
            n_total: total number of reads (ref + alt)
            genotype: the genotype of interest

        [Returns]
            the log-likelihood of observing n_ref, n_alt, given genotype
        """
        if genotype == 'R':
            result = log_betabinom_pmf(n_alt, n_total, self.alpha_R, self.beta_R)
        elif genotype == 'A':
            result = log_betabinom_pmf(n_alt, n_total, self.alpha_A, self.beta_A)
        elif genotype == 'H':
            result_no_dropout = np.log(1 - dropout_prob) + \
                                log_betabinom_pmf(n_alt, n_total, alpha_h, beta_h)
            dropout_R = np.log(dropout_prob) + np.log(1 - self.dropout_direction_prob) + \
                        log_betabinom_pmf(n_alt, n_total, self.alpha_R, self.beta_R)

            dropout_A = np.log(dropout_prob) + np.log(self.dropout_direction_prob) + \
                        log_betabinom_pmf(n_alt, n_total, self.alpha_A, self.beta_A)
            result = logsumexp([result_no_dropout, dropout_R, dropout_A])
        else:
            raise ValueError('[MutationFilter.single_read_llh] Invalid genotype.')

        return result


    def k_mut_llh(self, ref, alt, gt1, gt2):
        """
        [Arguments]
            ref, alt: 1D array, read counts at a locus for all cells
            gt1, gt2: genotypes before and after the mutation

        [Returns]
            If gt1 is the same as gt2 (i.e. there is no mutation), returns a single joint log-likelihood
            Otherwise, returns a 1D array in which entry [k] is the log-likelihood of having k mutated cells
        """

        N = ref.size  # number of cells
        total = ref + alt

        if gt1 == gt2:
            return np.sum([self.single_read_llh_with_dropout(alt[i], total[i], gt1) for i in range(N)])

        k_in_first_n_llh = np.zeros((N + 1, N + 1))  # [n,k]: log-likelihood that k among the first n cells are mutated
        k_in_first_n_llh[
            0, 0] = 0  # Trivial case: when there is 0 cell in total, the likelihood of having 0 mutated cell is 1

        for n in range(N):
            # log-likelihoods of the n-th cell having gt1 and gt2
            gt1_llh = self.single_read_llh_with_dropout(alt[n], total[n], gt1)
            gt2_llh = self.single_read_llh_with_dropout(alt[n], total[n], gt2)

            # k = 0 special case
            k_in_first_n_llh[n + 1, 0] = k_in_first_n_llh[n, 0] + gt1_llh

            # k = 1 through n
            k_over_n = np.array([k / (n + 1) for k in range(1, n + 1)])
            log_summand_1 = np.log(1 - k_over_n) + gt1_llh + k_in_first_n_llh[n, 1:n + 1]
            log_summand_2 = np.log(k_over_n) + gt2_llh + k_in_first_n_llh[n, 0:n]
            k_in_first_n_llh[n + 1, 1:n + 1] = np.logaddexp(log_summand_1, log_summand_2)

            # k = n+1 special case
            k_in_first_n_llh[n + 1, n + 1] = k_in_first_n_llh[n, n] + gt2_llh

        return k_in_first_n_llh[N, :]


    def single_locus_posteriors(self, ref, alt, comp_priors):
        """
        Calculates the log-posterior of different mutation types for a single locus

        # [Arguments]
        #     ref, alt: 1D arrays containing ref and alt reads of each cell
        #     comp_priors: log-prior for each genotype composition
        #
        # [Returns]
        #     1D numpy array containing posteriors of each mutation type, in the order
        #     ['R', 'H', 'A', 'RH', 'HA', 'HR', 'AH']
        #
        # NB When a mutation affects a single cell or all cells, it is considered non-mutated and assigned to one
        # of 'R', 'H' and 'A', depending on which one is the majority
        # """

        # f_R_est, f_H_est, f_A_est = self.em_algorithm_genotypes(ref, alt)
        # print(f_R_est, f_H_est, f_A_est)

        llh_RH = self.k_mut_llh(ref, alt, 'R', 'H')
        llh_HA = self.k_mut_llh(ref, alt, 'H', 'A')
        assert (llh_RH[-1] == llh_HA[0])  # both should be llh of all H

        joint_R = llh_RH[:1] + comp_priors[
            'R']  # llh zero out of n cells with genotype R are mutated given the data + prior of genotype RR
        joint_H = llh_HA[:1] + comp_priors['H']
        joint_A = llh_HA[-1:] + comp_priors['A']  # log likelihood, that all the cells are mutated + prior genotype A
        joint_RH = llh_RH[1:] + comp_priors[
            'RH']  # llh that 1 or more cells with genotype R are mutated given the data + prior of having 1 or more cells with genotype R having a mutation
        joint_HA = llh_HA[1:] + comp_priors['HA']
        joint_HR = np.flip(llh_RH)[1:] + comp_priors[
            'HR']  # llh k cells genotype H -> R  + prior that a mutation affects k cells for genotype H->R
        joint_AH = np.flip(llh_HA)[1:] + comp_priors['AH']

        joint = np.array([
            logsumexp(np.concatenate((joint_R, joint_RH[:1], joint_HR[-1:]))),  # RR
            logsumexp(np.concatenate((joint_H, joint_RH[-1:], joint_HR[:1], joint_HA[:1], joint_AH[-1:]))),  # HH
            logsumexp(np.concatenate((joint_A, joint_HA[-1:], joint_AH[:1]))),  # AA
            logsumexp(joint_RH[1:-1]),
            logsumexp(joint_HA[1:-1]),
            logsumexp(joint_HR[1:-1]),
            logsumexp(joint_AH[1:-1])
        ])

        posteriors = lognormalize(joint)  # Bayes' theorem

        return posteriors


    def mut_type_posteriors(self, ref, alt):
        """
        Calculates the log-prior of different mutation types for all loci
        In case no mutation occurs, all cells have the same genotype (which is either R or H or A)
        In case there is a mutation, each number of mutated cells is considered separately

        [Arguments]
            ref, alt: matrices containing the ref and alt reads

        [Returns]
            2D numpy array with n_loci rows and 7 columns, with each column standing for a mutation type
        """
        n_cells, n_loci = ref.shape

        # log-prior for each number of affected cells
        # placing a mutation randomly on one of the edges of a binary tree, how many affected cells would you expect?
        k_mut_priors = np.array([2 * logbinom(n_cells, k) - np.log(2 * k - 1) - logbinom(2 * n_cells, 2 * k) for k in
                                 range(1, n_cells + 1)])

        # composition priors
        comp_priors = {}
        for mut_type in ['R', 'H', 'A']:
            comp_priors[mut_type] = self.mut_type_prior[mut_type]
        for mut_type in ['RH', 'HA', 'HR', 'AH']:
            comp_priors[mut_type] = self.mut_type_prior[mut_type] + k_mut_priors

        # calculate posteriors for all loci with the help of multiprocessing
        result = np.zeros((n_loci, 7))
        for j in range(n_loci):
            result[j] = np.exp(self.single_locus_posteriors(ref[:, j], alt[:, j], comp_priors))
        return result


    def filter_mutations(self, ref, alt, method='highest_post', t=None, n_exp=None, ref_to_alt=False):
        """
        Filters the loci according to the posteriors of each mutation state

        [Arguments]
            method: criterion that determines which loci are considered mutated
            t: the posterior threshold to be used when using the 'threshold' method
            n_exp: the number of loci to be selected when using the 'first_k' method
        """
        assert (ref.shape == alt.shape)
        # ['R', 'H', 'A', 'RH', 'HA', 'HR', 'AH']
        posteriors = self.mut_type_posteriors(ref, alt)

        if method == 'highest_post':  # for each locus, choose the state with highest posterior
            selected = np.where(np.argmax(posteriors, axis=1) >= 3)[0]
        elif method == 'threshold':  # choose loci at which mutated posterior > threshold
            selected = np.where(np.sum(posteriors[:, 3:], axis=1) > t)[0]
        elif method == 'first_k':  # choose the k loci with highest mutated posteriors
            mut_posteriors = np.sum(posteriors[:, 3:], axis=1)
            order = np.argsort(mut_posteriors)[::-1]
            selected = order[:n_exp]
        else:
            raise ValueError('[MutationFilter.filter_mutations] Unknown filtering method.')

        if ref_to_alt: # If we only allow mutations that change a reference allele to an alternative allele
            mut_type = np.argmax(posteriors[selected, 3:5], axis=1)
            gt1_inferred = np.choose(mut_type, choices=['R', 'H'])
            gt2_inferred = np.choose(mut_type, choices=['H', 'A'])

        else:
            mut_type = np.argmax(posteriors[selected, 3:], axis=1)
            gt1_inferred = np.choose(mut_type, choices=['R', 'H', 'H', 'A'])
            gt2_inferred = np.choose(mut_type, choices=['H', 'A', 'R', 'H'])

        gt_not_selected = []  # maximum likelihood of genotypes that are not included in the tree learning
        for i in range(posteriors.shape[0]):
            if i in selected:
                continue
            else:
                genotype = np.argmax(posteriors[i, :3])
                gt_not_selected.append(np.choose(genotype, choices=['R', 'H', 'A']))

        return selected, gt1_inferred, gt2_inferred, gt_not_selected


    def get_llh_mat(self, ref, alt, gt1, gt2, individual=False, dropout_probs=None,
                    overdispersions_h=None):
        """
        [Arguments]
            individual: Use an individual dropout likelihood per SNV
        [Returns]
            llh_mat_1: 2D array in which [i,j] is the log-likelihood of cell i having gt1 at locus j
            llh_mat_2: 2D array in which [i,j] is the log-likelihood of cell i having gt2 at locus j
        """
        n_cells, n_mut = ref.shape
        llh_mat_1 = np.empty((n_cells, n_mut))
        llh_mat_2 = np.empty((n_cells, n_mut))
        total = ref + alt

        if individual:
            alphas_h = np.array(overdispersions_h) * 0.5
            betas_h = np.array(overdispersions_h) * 0.5
            assert n_mut == len(alphas_h) and n_mut == len(betas_h)

        for j in range(n_mut):
            for i in range(n_cells):
                k = alt[i, j]
                n = total[i, j]

                if not individual:
                    llh_mat_1[i, j] = self.single_read_llh_with_dropout(k, n, gt1[j]) # + self.mut_type_prior[gt1[j]]
                    llh_mat_2[i, j] = self.single_read_llh_with_dropout(k, n, gt2[j]) # + self.mut_type_prior[gt2[j]]
                else:
                    llh_mat_1[i, j] = self.single_read_llh_with_individual_dropout(k, n, gt1[j],
                                                                                   dropout_probs[j],
                                                                                   alphas_h[j], betas_h[j]) # + self.mut_type_prior[gt1[j]])
                    llh_mat_2[i, j] = self.single_read_llh_with_individual_dropout(k, n, gt2[j],
                                                                                   dropout_probs[j],
                                                                                   alphas_h[j], betas_h[j]) # + self.mut_type_prior[gt2[j]]

        return llh_mat_1, llh_mat_2


    def compute_llh_parameters(self, dropout_prob, overdispersion, error_rate, overdispersion_h,
                          global_opt=False, min_value=2, shape=2, alpha_parameters=2):
        """
        Computes the log-likelihood of the parameters of the model based on the priors.
        """

        alpha_dropout_prob = alpha_parameters
        beta_dropout_prob = 1 / self.dropout_prob

        alpha_error_rate = alpha_parameters
        beta_error_rate = 1 / self.error_rate

        divisor = shape # for individual parameter optimization use mean
        if global_opt: # for global parameter optimization use mode
            divisor = shape - 1
        scale_overdispersion_h = (self.overdispersion_H - min_value) / divisor
        scale_overdispersion = (self.overdispersion - min_value) / divisor

        return (
                beta.logpdf(dropout_prob, alpha_dropout_prob, beta_dropout_prob) +
                gamma.logpdf(overdispersion, shape, loc=min_value, scale=scale_overdispersion)+
                beta.logpdf(error_rate, alpha_error_rate, beta_error_rate) +
                gamma.logpdf(overdispersion_h, shape, loc=min_value, scale=scale_overdispersion_h)
        )


    def total_log_likelihood(self, params, k_obs, n_obs, genotypes):
        """
        Computes the total log-likelihood of the observations for the given parameters and genotypes.
        """
        dropout_prob, overdispersion, error_rate, overdispersion_h = params
        log_likelihood = 0

        # we assume, that the cells are independent in their allelic expression imbalance
        alpha_h = 0.5 * overdispersion_h
        beta_h = overdispersion_h - alpha_h

        for k, n, genotype in zip(k_obs, n_obs, genotypes):
            if genotype == "R":
                alpha_R = error_rate * overdispersion
                beta_R = overdispersion - alpha_R
                log_likelihood += log_betabinom_pmf(k, n, alpha_R, beta_R)

            elif genotype == "H":
                log_no_dropout, log_dropout_R, log_dropout_A = calculate_heterozygous_log_likelihoods(
                    k, n, dropout_prob, self.dropout_direction_prob, alpha_h, beta_h, error_rate, overdispersion)

                # Combine probabilities
                log_likelihood += logsumexp([log_no_dropout, log_dropout_R, log_dropout_A])

            elif genotype == "A":

                alpha_A = (1 - error_rate) * overdispersion
                beta_A = overdispersion - alpha_A
                log_likelihood += log_betabinom_pmf(k, n, alpha_A, beta_A)

            else:
                raise ValueError(f"Unexpected genotype: {genotype}")

        return log_likelihood


    def total_log_posterior(self, params, k_obs, n_obs, genotypes):
        """
        Computes the combined total log-likelihood of the parameters and observations.
        """
        dropout_prob, overdispersion, error_rate, overdispersion_h = params

        # Compute log-likelihood of the read counts given the parameters
        log_likelihood = self.total_log_likelihood(params, k_obs, n_obs, genotypes)

        # Compute log-likelihood of the parameters
        log_prior = self.compute_llh_parameters(dropout_prob, overdispersion, error_rate,
                                           overdispersion_h, global_opt=True)

        return -(log_likelihood + log_prior)  # Negative because we're minimizing


    def total_log_posterior_individual(self, params, k_obs, n_obs, overdispersion, error_rate, dropout_direction_prob):
        """
        Computes the combined total log-likelihood of the parameters and observations for individual SNVs.
        """
        dropout_prob, overdispersion_h = params

        # Compute log-likelihood of the read counts given the parameters
        log_likelihood = 0

        alpha_h = overdispersion_h * dropout_direction_prob
        beta_h = overdispersion_h - alpha_h

        for k, n in zip(k_obs, n_obs):
            log_no_dropout, log_dropout_R, log_dropout_A = calculate_heterozygous_log_likelihoods(
                k, n, dropout_prob, dropout_direction_prob, alpha_h, beta_h, error_rate, overdispersion)

            # Combine probabilities
            log_likelihood += logsumexp([log_no_dropout, log_dropout_R, log_dropout_A])

        # Compute the log-likelihood of the parameters
        log_prior = self.compute_llh_parameters(dropout_prob, overdispersion, error_rate,
                                           overdispersion_h, global_opt=False)

        return -(log_likelihood + log_prior)  # Negative because we're minimizing


    def fit_parameters(self, ref, alt, genotypes, initial_params=None,
                       max_iterations=50,
                       tolerance=1e-5):
        """
        Fits the dropout probability, overdispersion for homozygous genotypes, error rate, and overdispersion for
        heterozygous genotypes based on the observed read counts and genotypes.
        Uses a gradient-based optimization method to minimize the negative log-posterior.
        The parameters are bounded to ensure they remain within reasonable ranges.
        [Arguments]
            ref: 1D array of reference read counts
            alt: 1D array of alternative read counts
            genotypes: 1D array of inferred genotypes corresponding to the read counts
            initial_params: Initial guess for the parameters [dropout_prob, overdispersion, error rate, overdispersion_h]
            max_iterations: Maximum number of iterations for the optimization algorithm
            tolerance: Tolerance for convergence of the optimization algorithm
        [Returns]
            A 1D array of optimized parameters [dropout_prob, overdispersion, error_rate, overdispersion_h]
        """
        bounds = [
            (0.01, 0.99),  # dropout_prob heterozygous
            (2.5, 100),  # overdispersion for homozygous genotypes
            (0.001, 0.1),  # error_rate
            (2.5, 50),  # overdispersion_h for heterozygous genotypes
        ]

        total = alt + ref
        ind_nonzero = np.where(total != 0)[0]
        genotypes_nonzero = genotypes[ind_nonzero]
        alt_norm = alt[ind_nonzero]
        total_norm = total[ind_nonzero]

        def objective(params):
            return self.total_log_posterior(params, alt_norm, total_norm, genotypes_nonzero)

        result = minimize(
            objective,
            initial_params,
            method="L-BFGS-B",  # Gradient-based, bounded optimization
            bounds=bounds,
            options={"maxiter": max_iterations, "ftol": tolerance}
        )

        if not result.success:
            print(f"Optimization failed: {result.message}")

        return result.x


    def fit_parameters_individual(self, alt_het, total_reads, overdispersion, error_rates, dropout_direction,
                                  initial_params=None, max_iterations=50, tolerance=1e-5):
        """
        Fits the dropout probability and overdispersion for heterozygous genotypes for individual SNVs.
        [Arguments]
            alt_het: 1D array of alternative read counts for heterozygous genotypes
            total_reads: 1D array of total read counts (ref + alt) for heterozygous genotypes
            overdispersion: Overdispersion for homozygous genotypes
            error_rates: Expected alternative allele frequency for homozygous reference genotypes
            dropout_direction: Frequency, that the reference allele is dropped out instead of the alternative in the
            heterozygous case
            initial_params: Initial guess for the parameters [dropout_prob, overdispersion_h]
            max_iterations: Maximum number of iterations for the optimization algorithm
            tolerance: Tolerance for convergence of the optimization algorithm
        [Returns]
            A 1D array of optimized parameters [dropout_prob, overdispersion_h]
        """
        bounds = [
            (0.01, 0.99),  # dropout_prob heterozygous
            (2.5, 50),  # overdispersion_h for heterozygous
        ]

        def objective(params):
            return self.total_log_posterior_individual(params, alt_het, total_reads, overdispersion, error_rates,
                                                       dropout_direction)

        result = minimize(
            objective,
            initial_params,
            method="L-BFGS-B",  # Gradient-based, bounded optimization
            bounds=bounds,
            options={"maxiter": max_iterations, "ftol": tolerance}
        )

        if not result.success:
            print(f"Optimization failed: {result.message}")
        return result.x


    def update_parameters(self, ref, alt, inferred_genotypes):
        """
        Update the parameters based on the inferred genotypes and observed read counts.
        Globally optimized parameters over all SNVs are the dropout probability,
        overdispersion (heterozygous genotypes), overdispersion (homozygous genotypes) and error rate. Parameters
        additionally estimated for each SNV with sufficient data are the dropout probability and
        overdispersion (heterozygous genotypes).
        """
        all_ref_counts = ref.flatten()
        all_alt_counts = alt.flatten()
        all_genotypes = inferred_genotypes.flatten()

        global_params = self.fit_parameters(all_ref_counts, all_alt_counts, all_genotypes,
                                            [self.dropout_prob, self.overdispersion, self.error_rate, self.overdispersion_H])
        print(f"Global parameters: {global_params}")

        dropout_prob, overdispersion, error_rate, overdispersion_h = global_params

        individual_dropout_probs = []
        individual_overdispersions_h = []

        for snv in range(len(inferred_genotypes[0])):
            indices = np.where(inferred_genotypes[:, snv] == "H")[0]

            ref_het = ref[indices, snv]
            alt_het = alt[indices, snv]

            total_reads = ref_het + alt_het

            informative = total_reads > 10  # Only use SNVs with more than 10 reads
            alt_het = alt_het[informative]
            total_reads = total_reads[informative]

            if len(total_reads) > 5:  # Fit individual parameters if there are enough informative heterozygous cells
                individual_params = self.fit_parameters_individual(
                    alt_het, total_reads, overdispersion, error_rate,
                    dropout_direction=self.dropout_direction_prob,
                    initial_params=[dropout_prob, overdispersion_h],
                )
            else:
                individual_params = dropout_prob, overdispersion_h

            posterior_dropout_prob, posterior_overdispersion_hs = individual_params

            individual_dropout_probs.append(posterior_dropout_prob)
            individual_overdispersions_h.append(posterior_overdispersion_hs)

        return dropout_prob, overdispersion, error_rate, overdispersion_h, \
            individual_dropout_probs, individual_overdispersions_h


# Additional helper functions
def calculate_heterozygous_log_likelihoods(k, n, dropout_prob, dropout_direction_prob, alpha_h, beta_h,
                                           error_rate, overdispersion):
    """
    Calculates the log-likelihood of observing k alternative reads with n coverage for a heterozygous locus.
    """
    log_no_dropout = np.log(1 - dropout_prob) + log_betabinom_pmf(k, n, alpha_h, beta_h)

    # Dropout to "R"
    alpha_R = error_rate * overdispersion
    beta_R = overdispersion - alpha_R
    log_dropout_R = np.log(dropout_prob) + np.log(1 - dropout_direction_prob) + log_betabinom_pmf(k, n, alpha_R, beta_R)

    # Dropout to "A"
    alpha_A = (1 - error_rate) * overdispersion
    beta_A = overdispersion - alpha_A
    log_dropout_A = np.log(dropout_prob) + np.log(dropout_direction_prob) + log_betabinom_pmf(k, n, alpha_A, beta_A)

    return log_no_dropout, log_dropout_R, log_dropout_A

@njit
def ln_gamma(z):
    return math.lgamma(z)


@njit
def betaln(x, y):
    return ln_gamma(x) + ln_gamma(y) - ln_gamma(x + y)


@njit
def factorial(n):
    if n == 0 or n == 1:
        return 0.0  # log(1) = 0, log(0) is undefined, return 0.0 as safe placeholder
    result = 0.0
    for i in range(2, n + 1):
        result += np.log(i)
    return result


def logbinom(n, k):
    return loggamma(n + 1) - loggamma(k + 1) - loggamma(n - k + 1)


def lognormalize(array):
    return array - logsumexp(array)


@njit
def log_binomial_coefficient(n, k):
    if 0 <= k <= n:
        log_numerator = factorial(n)
        log_denominator = factorial(k) + factorial(n - k)
        return log_numerator - log_denominator
    else:
        return 0.0  # Return 0.0 for invalid inputs

@njit
def log_betabinom_pmf(k, n, a, b):
    if n < 0 or k < 0 or k > n or a <= 0 or b <= 0:
        return 0.0

    log_binom_coef = log_binomial_coefficient(n, k)

    num = betaln(k + a, n - k + b)
    denom = betaln(a, b)

    return num - denom + log_binom_coef
