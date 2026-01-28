"""
This script is used to generate simulated read counts, ground truth genotypes and ground truth trees.
"""

import numpy as np
from scipy.stats import poisson, geom, nbinom, gamma

from src_python.cell_tree import CellTree
from src_python.mutation_tree import MutationTree
from src_python.utils import load_config_and_set_random_seed

config = load_config_and_set_random_seed()


def betabinom_rvs(coverage, alpha, beta_param):
    """
    Generate random samples from a Beta-Binomial distribution.
    """
    p = np.random.beta(alpha, beta_param)
    return np.random.binomial(coverage, p)


class DataGenerator:
    """
    DataGenerator is a class used to generate simulated read counts, ground truth genotypes, and ground truth trees
    for single-cell sequencing data. Mutation types are ["RH", "AH", "HR", "HA"]

    Attributes:
        n_cells – Number of cells.
        n_mut – Number of mutations.
        mut_prop – Proportion of mutations.
        error_rate – Expected variant allele frequency for homozygous reference genotype
        overdispersion – Overdispersion parameter in homozygous reference and alternative cases
        genotype_freq – Frequencies of genotypes.
        coverage_method – Method to generate coverage values.
        coverage_mean – Mean coverage value.
        coverage_sample – Array of coverage values to sample from.
        dropout_alpha – Alpha parameter of the beta distribution the dropout probability is sampled from.
        dropout_beta – Beta parameter of the beta distribution the dropout probability is sampled from.
        dropout_dir – Dropout direction probability.
        overdispersion_h – Overdispersion parameter for the read counts in the heterozygous case.
        homoplasy_fraction - Fraction of loci with two mutations per locus in different branches of the tree
        CNV_fraction - Fraction of loci affected by copy number alterations
    """

    def __init__(self, n_cells, n_mut, coverage_method="zinb", genotype_freq=None, coverage_sample=None,
                 dropout_alpha=None, dropout_beta=None, **kwargs):

        self.ref_alleles = None
        self.alt_alleles = None
        self.dropout = kwargs.get("dropout", 0.2)
        self.overdispersion_h = kwargs.get("overdispersion_Het", 6)
        self.overdispersion = kwargs.get("overdispersion_Hom", 10)
        self.error_rate = kwargs.get("error_rate", 0.05)
        self.coverage_mean = kwargs.get("coverage_mean", 60)
        self.coverage_zero_inflation = kwargs.get("coverage_zero_inflation", 0.39)
        self.coverage_dispersion = kwargs.get("coverage_dispersion", 5.88)
        self.homoplasy_fraction = kwargs.get("homoplasy_fraction", 0.0)
        self.CNV_fraction = kwargs.get("CNV_fraction", 0.0)

        if dropout_alpha is not None and dropout_beta is not None:
            self.dropout_alpha = dropout_alpha
            self.dropout_beta = dropout_beta
        else:
            self.dropout_alpha = self.dropout * 10
            self.dropout_beta = 10 - self.dropout_alpha

        self.dropout_dir = config["dropout_direction"]

        self.ct = CellTree(n_cells=n_cells, n_mut=n_mut)
        self.mt = MutationTree(n_mut=n_mut, n_cells=n_cells)

        self.genotype = np.empty((self.n_cells, self.n_mut), dtype=str)
        self.mut_prop = 1
        self.genotype_freq = [1 / 3, 1 / 3, 1 / 3] if genotype_freq is None else genotype_freq
        self.gt1 = np.random.choice(["R", "H", "A"], size=self.n_mut, replace=True, p=self.genotype_freq)
        self.gt2 = np.empty_like(self.gt1)

        self.coverage_method = coverage_method

        if coverage_method == "sample" and coverage_sample is None:
            raise ValueError("Please provide array of coverage values to be sampled from.")
        self.coverage_sample = coverage_sample

        self.alpha_R = self.error_rate * self.overdispersion
        self.beta_R = self.overdispersion - self.alpha_R
        self.alpha_A = (1 - self.error_rate) * self.overdispersion
        self.beta_A = self.overdispersion - self.alpha_A

    @property
    def n_cells(self):
        return self.ct.n_cells

    @property
    def n_mut(self):
        return self.ct.n_mut

    def random_mut_type(self):
        self.gt1 = np.random.choice(["R", "H", "A"], size=self.n_mut, replace=True, p=self.genotype_freq)
        mutated = np.random.choice(self.n_mut, size=round(self.n_mut * self.mut_prop), replace=False)
        for j in mutated:
            if self.gt1[j] == "H":
                self.gt2[j] = np.random.choice(["R", "A"])  # mutation HA and HR with equal probability
            else:
                self.gt2[j] = "H"

    def random_tree(self, num_clones):
        """
        Generate a random tree

        [Arguments]
            num_clones: Number of clones (cells with the same genotype) in the tree.
        """
        self.ct.rand_subtree()
        self.ct.rand_mut_loc(num_clones)

    def random_coverage(self):
        """
        Generate random coverage values for each cell and mutation based on the specified coverage method.

        The coverage values are generated using one of the following methods:
        - "constant": All coverage values are set to the mean coverage.
        - "geometric": Coverage values are sampled from a geometric distribution.
        - "poisson": Coverage values are sampled from a Poisson distribution.
        - "zinb": Coverage values are sampled from a zero-inflated negative binomial distribution.
        - "sample": Coverage values are sampled from a provided array of coverage values.
        """
        match self.coverage_method:
            case "constant":
                self.coverage = np.ones((self.n_cells, self.n_mut), dtype=int) * self.coverage_mean
            case "geometric":
                self.coverage = geom.rvs(p=1 / (self.coverage_mean + 1), loc=-1, size=(self.n_cells, self.n_mut))
            case "poisson":
                self.coverage = poisson.rvs(mu=self.coverage_mean, size=(self.n_cells, self.n_mut))
            # parameters 60, 5.88, 0.38 learned from mm34 scRNA seq dataset
            case "zinb":
                mu, alpha, pi = self.coverage_mean, self.coverage_dispersion, self.coverage_zero_inflation # alpha: overdispersion, pi: zero-inflation probability
                n = 1 / alpha  # number of successes
                p = 1 / (1 + alpha * mu)
                nb_samples = nbinom.rvs(n, p, size=(self.n_cells, self.n_mut))
                zero_inflation_mask = np.random.rand(self.n_cells, self.n_mut) < pi
                self.coverage = np.where(zero_inflation_mask, 0, nb_samples)
            case "sample":
                self.coverage = np.random.choice(self.coverage_sample, size=(self.n_cells, self.n_mut), replace=True)
            case _:
                raise ValueError("Invalid coverage sampling method.")

    def generate_single_read(self, ref_alleles, alt_alleles, genotype, coverage, dropout_prob, dropout_direction, overdispersion_h):

        if genotype == "R":
            n_alt = betabinom_rvs(coverage, self.alpha_R, self.beta_R)
        elif genotype == "A":
            n_alt = betabinom_rvs(coverage, self.alpha_A, self.beta_A)
        elif genotype == "H":
            # Determine if dropout occurs
            dropout_occurs = np.random.rand() < dropout_prob

            if dropout_occurs:
                if np.random.rand() < dropout_direction:
                    ref_alleles -= 1
                else:
                    alt_alleles -= 1

            if ref_alleles == 0:
                # Only alt alleles remain; dropout to A
                n_alt = betabinom_rvs(coverage, self.alpha_A, self.beta_A)
            elif alt_alleles == 0:
                # Only ref alleles remain; dropout to R
                n_alt = betabinom_rvs(coverage, self.alpha_R, self.beta_R)
            else:
                # Both alleles still present
                cna = (ref_alleles + alt_alleles)
                alpha_h = (alt_alleles / cna) * overdispersion_h * cna # scale with copy number to maintain hill shape
                beta_h = overdispersion_h * cna - alpha_h
                n_alt = betabinom_rvs(coverage, alpha_h, beta_h)
        else:
            raise ValueError("[generate_single_read] ERROR: invalid genotype.")

        n_ref = coverage - n_alt
        return n_ref, n_alt


    def generate_reads(self, new_tree=False, new_mut_type=False, new_coverage=True, num_clones="", min_value=2.5,
                       shape=2):
        """
        Generate read counts for all cells and mutations.
        [Arguments]
            new_tree: if True, generate a new random tree
            new_mut_type: if True, generate new random mutation types
            new_coverage: if True, generate new random coverage values
            num_clones: number of clones in the tree, otherwise random mutation placement is used
            min_value: minimum value for the overdispersion parameter in the heterozygous case
            shape: shape parameter for the gamma distribution, which is used to sample the overdispersion parameter

        [Returns]
            ref: reference read counts
            alt: alternative read counts
            all_dropout_probs: dropout probabilities per SNV
            all_overdispersions_h: overdispersion parameters for heterozygous case per SNV
        """
        if new_mut_type:
            self.random_mut_type()
        if new_tree:
            self.random_tree(num_clones)
        if new_coverage:
            self.random_coverage()

        # track the number of ref and alt alleles for CNV simulation
        self.ref_alleles = np.zeros((self.n_cells, self.n_mut), dtype=int)
        self.alt_alleles = np.zeros((self.n_cells, self.n_mut), dtype=int)

        # determine genotypes
        self.genotype = np.empty((self.n_cells, self.n_mut), dtype=str)
        mut_indicator = self.mut_indicator()
        for i in range(self.n_cells):
            for j in range(self.n_mut):
                self.genotype[i, j] = self.gt2[j] if mut_indicator[i, j] else self.gt1[j]

                if self.genotype[i, j] == "R":
                    self.ref_alleles[i, j], self.alt_alleles[i, j] = 2, 0
                elif self.genotype[i, j] == "H":
                    self.ref_alleles[i, j], self.alt_alleles[i, j] = 1, 1
                elif self.genotype[i, j] == "A":
                    self.ref_alleles[i, j], self.alt_alleles[i, j] = 0, 2

        # Apply CNVs
        for j in range(self.n_mut):
            for i in range(self.n_cells):
                if np.random.random() < self.CNV_fraction:
                    cnv = np.random.choice([1, 3, 4, 5, 6])

                    self.coverage[i, j] = int(cnv/2 * self.coverage[i, j]) # adjust coverage according to CNV

                    current_alleles = []
                    if self.genotype[i, j] == "R":
                        current_alleles = ["ref", "ref"]
                    elif self.genotype[i, j] == "H":
                        current_alleles = ["ref", "alt"]
                    elif self.genotype[i, j] == "A":
                        current_alleles = ["alt", "alt"]

                    if cnv == 1:
                        # Randomly drop one allele
                        if current_alleles:
                            dropped_allele = np.random.choice(current_alleles)
                            if dropped_allele == "ref":
                                self.ref_alleles[i, j] = self.ref_alleles[i, j] - 1
                            else:
                                self.alt_alleles[i, j] = self.alt_alleles[i, j] - 1
                    else:
                        # For CNVs > 1, duplicate alleles (cnv-2) times
                        for _ in range(cnv - 2):
                            chosen_allele = np.random.choice(current_alleles)
                            if chosen_allele == "ref":
                                self.ref_alleles[i, j] += 1
                            else:
                                self.alt_alleles[i, j] += 1
                            # Update current_alleles to reflect the new allele count
                            current_alleles.append(chosen_allele)

        for i in range(self.n_cells):
            for j in range(self.n_mut):
                if self.alt_alleles[i, j] == 0:
                    self.genotype[i, j] = "R"
                elif self.ref_alleles[i, j] == 0:
                    self.genotype[i, j] = "A"
                else:
                    self.genotype[i, j] = "H"


        # read count generation
        ref = np.empty((self.n_cells, self.n_mut), dtype=int)
        alt = np.empty((self.n_cells, self.n_mut), dtype=int)

        all_dropout_probs = []
        all_overdispersions_h = []

        for j in range(self.n_mut):
            # Sample dropout probabilities from beta distributions for each SNV
            if self.dropout != 0:
                dropout_prob = np.random.beta(self.dropout_alpha, self.dropout_beta)
            else:
                dropout_prob = 0.0

            all_dropout_probs.append(dropout_prob)

            # Sample overdispersion parameter for heterozygous case from gamma distribution for each SNV
            scale = (self.overdispersion_h - min_value) / shape  # θ = (mean - shift) / k

            overdispersion_H = gamma.rvs(shape, loc=min_value, scale=scale)
            all_overdispersions_h.append(overdispersion_H)

            # as the cells are assumed to be independent, allelic imbalances
            # are assumed to be symmetric and only affect the overdispersion

            for i in range(self.n_cells):
                ref[i, j], alt[i, j] = self.generate_single_read(self.ref_alleles[i, j], self.alt_alleles[i,j], self.genotype[i, j], self.coverage[i, j],
                                                                 dropout_prob, config["dropout_direction"], overdispersion_H)

        return ref, alt, all_dropout_probs, all_overdispersions_h

    def mut_indicator(self):
        """
        Return a 2D Boolean array in which [i,j] indicates whether cell i is affected by mutation j
        """
        res = np.zeros((self.n_cells, self.n_mut), dtype=bool)
        for j in range(self.n_mut):  # determine for each mutation the cells below it in the tree
            for i in self.ct.leaves(self.ct.mut_loc[j]):
                res[i, j] = True

            if np.random.rand() < self.homoplasy_fraction:
                # place the mutation a second time independently (homoplasy)
                loc1 = self.ct.mut_loc[j]
                ancestors = [a for a in self.ct.ancestors(loc1)]
                descendants = [d for d in self.ct.dfs(loc1)]
                options = [o for o in range(len(self.ct.parent_vec)) if o not in ancestors and o not in descendants]

                if len(options) == 0:
                    continue  # no valid location for the second placement

                loc2 = np.random.choice(options)

                for i in self.ct.leaves(loc2):
                    res[i, j] = True
        return res
