use itertools::Itertools;
use rand::{seq::IteratorRandom, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;

use crate::{
    amalgam::Amalgam,
    config::{TransformConfig, WeightConfig},
    math::{aggregate_pvalues, arithmetic_mean, empirical_fdr},
    results::{GeneResult, GeoPAGGResults},
    utils::{calculate_group_sizes, index_mask, select_indices},
};

/// Implementation of the GeoPAGG (Geometric P-Value Aggregation for Gene Grouping) Algorithm
pub struct GeoPAGG<'a> {
    pvalues: Vec<f64>,
    logfc: &'a [f64],
    genes: &'a [String],
    token: Option<&'a str>,
    weight_config: WeightConfig,
    seed: usize,
}

impl<'a> GeoPAGG<'a> {
    /// Creates a new GeoPAGG instance
    ///
    /// # Arguments
    ///
    /// * `pvalues` - Slice of p-values for each sgRNA
    /// * `logfc` - Slice of log fold changes for each sgRNA
    /// * `genes` - Slice of gene names for each sgRNA
    /// * `token` - Optional token for filtering genes
    /// * `weight_config` - Configuration for p-value weighting
    /// * `transform_config` - Configuration for p-value transformation
    /// * `seed` - Seed for random number generation
    #[must_use]
    pub fn new(
        pvalues: &'a [f64],
        logfc: &'a [f64],
        genes: &'a [String],
        token: Option<&'a str>,
        weight_config: WeightConfig,
        transform_config: TransformConfig,
        seed: usize,
    ) -> Self {
        let pvalues = transform_config.transform(pvalues);
        Self {
            pvalues,
            logfc,
            genes,
            token,
            weight_config,
            seed,
        }
    }

    /// Runs the GeoPAGG algorithm
    ///
    /// This method performs the following steps:
    /// 1. Build the amalgam groups
    /// 2. Aggregate each gene
    /// 3. Aggregate the amalgam groups
    /// 4. Combine the results and calculate the empirical FDR
    ///
    /// # Returns
    ///
    /// A `GeoPAGGResults` struct containing the analysis results
    pub fn run(&self) -> GeoPAGGResults {
        let unique_genes = self.genes.iter().unique().collect::<Vec<_>>();

        // Build the amalgam groups
        let mut rng = ChaCha8Rng::seed_from_u64(self.seed as u64);
        let null_set = self.distinguish_null_set();
        let group_sizes = calculate_group_sizes(self.genes, &unique_genes);
        let amalgams = self.build_amalgams(&group_sizes, &null_set, &mut rng);

        // Aggregate each gene
        let test_results = unique_genes
            .par_iter()
            .map(|gene| self.process_gene(gene))
            .collect::<Vec<_>>();

        // Aggregate the amalgams
        let amalgam_results = amalgams
            .into_par_iter()
            .map(|amalgam| amalgam.into())
            .collect::<Vec<_>>();

        // Combine the results
        let mut results = test_results
            .into_iter()
            .chain(amalgam_results)
            .collect::<Vec<_>>();

        empirical_fdr(&mut results);

        GeoPAGGResults::from_vec(results)
    }

    /// Returns the indices of all sgRNAs that can be incorporated into the amalgam groups
    ///
    /// If a token is provided, only sgRNAs whose genes match the token are included.
    /// Otherwise, all sgRNAs are considered for amalgam inclusion.
    fn distinguish_null_set(&self) -> Vec<usize> {
        let mut null_set = Vec::new();
        if let Some(token) = &self.token {
            let token_indices = index_mask(token, self.genes);
            null_set.extend(token_indices);
        } else {
            null_set.extend(0..self.genes.len());
        }

        null_set
    }

    /// Builds the amalgam groups
    ///
    /// Amalgam groups are built by randomly selecting sgRNAs from the null set,
    /// matching the membership sizes and distributions of the unique genes.
    fn build_amalgams(
        &self,
        group_sizes: &Vec<(usize, usize)>,
        null_set: &[usize],
        rng: &mut ChaCha8Rng,
    ) -> Vec<Amalgam> {
        let mut amalgams = Vec::new();
        for (membership_size, num_amalgams) in group_sizes {
            for draw_index in 0..*num_amalgams {
                let random_indices = null_set
                    .iter()
                    .copied()
                    .choose_multiple(rng, *membership_size);
                let pvalues = select_indices(&random_indices, &self.pvalues);
                let logfc = select_indices(&random_indices, self.logfc);
                let amalgam = Amalgam::new(
                    *membership_size,
                    draw_index,
                    pvalues,
                    logfc,
                    self.weight_config,
                );
                amalgams.push(amalgam);
            }
        }
        amalgams
    }

    /// Processes a single gene
    ///
    /// This function aggregates the p-values and log fold changes for a single gene
    /// and returns the results in a `GeneResult` struct.
    fn process_gene(&self, gene: &str) -> GeneResult {
        let gene_indices = index_mask(gene, self.genes);
        let mut pvalues = select_indices(&gene_indices, &self.pvalues);
        let logfc = select_indices(&gene_indices, self.logfc);

        // Aggregate the p-values
        let wgm = aggregate_pvalues(&mut pvalues, self.weight_config);

        // Aggregate the log fold changes as standard mean
        let logfc = arithmetic_mean(&logfc);

        GeneResult::new(gene.to_string(), wgm, logfc, false)
    }
}
