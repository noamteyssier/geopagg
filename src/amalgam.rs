use derive_new::new;

use crate::{
    config::WeightConfig,
    math::{aggregate_pvalues, arithmetic_mean},
    results::GeneResult,
};

/// Represents an amalgam group in the GeoPAGG algorithm
///
/// An amalgam group is a collection of randomly selected sgRNAs used for
/// comparison with actual gene groups to calculate empirical FDR.
#[derive(new)]
pub struct Amalgam {
    /// Number of sgRNAs in this amalgam group
    membership_size: usize,
    /// Index of this amalgam within its size group
    draw_index: usize,
    /// P-values of the sgRNAs in this amalgam
    pvalues: Vec<f64>,
    /// Log fold changes of the sgRNAs in this amalgam
    logfc: Vec<f64>,
    /// Configuration for p-value weighting
    weight_config: WeightConfig,
}

impl From<Amalgam> for GeneResult {
    /// Converts an Amalgam into a GeneResult
    ///
    /// This conversion aggregates the p-values and log fold changes of the amalgam
    /// and creates a unique identifier for the amalgam.
    fn from(amalgam: Amalgam) -> Self {
        let mut pvalues = amalgam.pvalues;
        let wgm = aggregate_pvalues(&mut pvalues, amalgam.weight_config);
        let logfc = arithmetic_mean(&amalgam.logfc);
        let gene = format!("amalgam_{}_{}", amalgam.membership_size, amalgam.draw_index);
        GeneResult::builder()
            .gene(gene)
            .wgm(wgm)
            .logfc(logfc)
            .amalgam(true)
            .build()
    }
}
