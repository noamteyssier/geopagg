use derive_new::new;

use crate::{
    config::WeightConfig,
    math::{aggregate_pvalues, arithmetic_mean},
    results::GeneResult,
};

/// An amalgam group
#[derive(new)]
pub struct Amalgam {
    membership_size: usize,
    draw_index: usize,
    pvalues: Vec<f64>,
    logfc: Vec<f64>,
    weight_config: WeightConfig,
}

impl From<Amalgam> for GeneResult {
    fn from(amalgam: Amalgam) -> Self {
        let mut pvalues = amalgam.pvalues;
        let wgm = aggregate_pvalues(&mut pvalues, amalgam.weight_config);
        let logfc = arithmetic_mean(&amalgam.logfc);
        let gene = format!("amalgam_{}_{}", amalgam.membership_size, amalgam.draw_index);
        GeneResult {
            gene,
            wgm,
            logfc,
            amalgam: true,
            empirical_fdr: 1.0,
            adjusted_empirical_fdr: 1.0,
        }
    }
}
