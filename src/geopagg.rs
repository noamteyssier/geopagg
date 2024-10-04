use bon::bon;
use itertools::Itertools;
use rand::{seq::IteratorRandom, SeedableRng};
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;

use crate::{
    amalgam::Amalgam,
    config::{TransformConfig, WeightConfig},
    math::{aggregate_pvalues, arithmetic_mean, empirical_fdr, empirical_fdr_product},
    results::{GeneResult, GeoPAGGResults},
    utils::{
        calculate_group_sizes, index_mask, index_prefix_mask, select_indices, zscore_transform,
    },
};

/// Implementation of the GeoPAGG (Geometric P-Value Aggregation for Gene Grouping) Algorithm
pub struct GeoPAGG<'a> {
    pvalues: Vec<f64>,
    logfc: &'a [f64],
    genes: &'a [String],
    token: Option<&'a str>,
    weight_config: WeightConfig,
    use_product: bool,
    ascending: bool,
    seed: usize,
    zscore_threshold: Option<f64>,
}

#[bon]
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
    #[builder]
    pub fn new(
        pvalues: &'a [f64],
        logfc: &'a [f64],
        genes: &'a [String],
        token: Option<&'a str>,
        #[builder(default)] //
        weight_config: WeightConfig,
        #[builder(default)] //
        transform_config: TransformConfig,
        #[builder(default)] //
        use_product: bool,
        #[builder(default)] //
        ascending: bool,
        seed: usize,
        zscore_threshold: Option<f64>,
    ) -> Self {
        let pvalues = transform_config.transform(pvalues);
        Self {
            pvalues,
            logfc,
            genes,
            token,
            weight_config,
            use_product,
            ascending,
            seed,
            zscore_threshold,
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
        let mut unique_genes = self.genes.iter().unique().collect::<Vec<_>>();
        unique_genes.sort();

        // Build the amalgam groups
        let mut rng = ChaCha8Rng::seed_from_u64(self.seed as u64);
        let mut null_set = self.distinguish_null_set();
        if let Some(zscore_threshold) = self.zscore_threshold {
            null_set = self.filter_null_set(&null_set, &self.pvalues, zscore_threshold);
        }
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

        // Calculate the empirical FDR
        if self.use_product {
            // Using the product of p-values and logfc can be thought of as stepping diagonally on the volcano plot.
            empirical_fdr_product(&mut results, self.ascending)
        } else {
            // Using the p-values by default can be thought of as stepping vertically on the volcano plot.
            empirical_fdr(&mut results)
        }

        GeoPAGGResults::from_vec(results)
    }

    /// Returns the indices of all sgRNAs that can be incorporated into the amalgam groups
    ///
    /// If a token is provided, only sgRNAs whose genes match the token are included.
    /// Otherwise, all sgRNAs are considered for amalgam inclusion.
    fn distinguish_null_set(&self) -> Vec<usize> {
        let mut null_set = Vec::new();
        if let Some(token) = &self.token {
            let token_indices = index_prefix_mask(token, self.genes);
            null_set.extend(token_indices);
        } else {
            null_set.extend(0..self.genes.len());
        }

        null_set
    }

    /// Filters the null set based on the z-score threshold
    fn filter_null_set(
        &self,
        null_set: &[usize],
        pvalues: &[f64],
        zscore_threshold: f64,
    ) -> Vec<usize> {
        let null_pvalues = select_indices(null_set, pvalues);
        let neg_log_pvalues: Vec<_> = null_pvalues.iter().map(|x| x.ln()).map(|x| -x).collect();
        let zscores = zscore_transform(&neg_log_pvalues);
        null_set
            .iter()
            .zip(zscores)
            .filter(|(_, zscore)| zscore.abs() < zscore_threshold)
            .map(|(index, _)| *index)
            .collect()
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

        GeneResult::builder()
            .gene(gene.to_string())
            .wgm(wgm)
            .logfc(logfc)
            .build()
    }
}

#[cfg(test)]
mod testing {
    use super::*;

    #[test]
    fn test_geopagg() {
        let pvalues = vec![0.1, 0.2, 0.3, 0.4, 0.5];
        let logfc = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let genes = vec![
            "A".to_string(),
            "B".to_string(),
            "A".to_string(),
            "B".to_string(),
            "A".to_string(),
        ];
        let seed = 42;

        let geopagg = GeoPAGG::builder()
            .pvalues(&pvalues)
            .logfc(&logfc)
            .genes(&genes)
            .seed(seed)
            .build();

        let results = geopagg.run();

        // 2 unique genes, 2 amalgam groups
        assert_eq!(results.empirical_fdr.len(), 4);

        // Checks that the results are sorted by gene name
        for i in 0..results.genes.len() - 1 {
            assert!(results.genes[i] <= results.genes[i + 1]);
        }

        // Checks that the proper genes are present
        // A, B, amalgam_3_0, amalgam_2_0
        for gene in &["A", "B", "amalgam_3_0", "amalgam_2_0"] {
            assert!(results.genes.contains(&gene.to_string()));
        }
        assert_eq!(results.genes.len(), 4);
    }

    #[test]
    fn test_seeding() {
        let pvalues = vec![0.1, 0.2, 0.3, 0.4, 0.5];
        let logfc = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let genes = vec![
            "A".to_string(),
            "B".to_string(),
            "A".to_string(),
            "B".to_string(),
            "A".to_string(),
        ];

        // Seed 42
        // Should be deterministic
        let mut last_42 = vec![];
        for idx in 0..1000 {
            let seed = 42;
            let geopagg = GeoPAGG::builder()
                .pvalues(&pvalues)
                .logfc(&logfc)
                .genes(&genes)
                .seed(seed)
                .build();
            let results = geopagg.run();
            if idx > 0 {
                assert_eq!(results.adjusted_empirical_fdr, last_42);
            }
            last_42 = results.adjusted_empirical_fdr;
        }

        // Seed 0
        // also deterministic
        let mut last_0 = vec![];
        for idx in 0..1000 {
            let seed = 0;
            let geopagg = GeoPAGG::builder()
                .pvalues(&pvalues)
                .logfc(&logfc)
                .genes(&genes)
                .seed(seed)
                .build();
            let results = geopagg.run();
            if idx > 0 {
                dbg!(idx);
                dbg!(&last_0);
                dbg!(&results.adjusted_empirical_fdr);
                assert_eq!(results.adjusted_empirical_fdr, last_0);
            }
            last_0 = results.adjusted_empirical_fdr;
        }
        // Test that the results are deterministic
        for seed in 0..100 {
            let mut last_seed = vec![];
            for idx in 0..50 {
                let geopagg = GeoPAGG::builder()
                    .pvalues(&pvalues)
                    .logfc(&logfc)
                    .genes(&genes)
                    .seed(seed)
                    .build();
                let results = geopagg.run();
                if idx > 0 {
                    assert_eq!(results.adjusted_empirical_fdr, last_seed);
                }
                last_seed = results.adjusted_empirical_fdr;
            }
        }
    }
}
