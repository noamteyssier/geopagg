use crate::{config::WeightConfig, results::GeneResult};

pub fn weighted_geometric_mean(x: &[f64], weights: &[f64]) -> f64 {
    let sum = x
        .iter()
        .zip(weights.iter())
        .map(|(x, w)| x.ln() * w)
        .sum::<f64>();
    let sum_of_weights = weights.iter().sum::<f64>();
    (sum / sum_of_weights).exp()
}

pub fn arithmetic_mean(x: &[f64]) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}

pub fn aggregate_pvalues(pvalues: &mut [f64], weight_config: WeightConfig) -> f64 {
    // Sort the pvalues
    pvalues.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

    let weights = weight_config.build_weights(pvalues.len());
    weighted_geometric_mean(pvalues, &weights)
}

/// Calculates the empirical FDR for each gene
/// by counting the number of amalgams that are lower than the gene's WGM
pub fn empirical_fdr(gene_results: &mut [GeneResult]) {
    gene_results.sort_unstable_by(|a, b| a.wgm.partial_cmp(&b.wgm).unwrap());

    let mut amalgam_count = 0;
    for (rank, gene_result) in gene_results.iter_mut().enumerate() {
        if gene_result.amalgam {
            amalgam_count += 1;
        }
        gene_result.empirical_fdr = amalgam_count as f64 / (rank + 1) as f64;
        gene_result.adjusted_empirical_fdr = gene_result.empirical_fdr.max(gene_result.wgm);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_weighted_geometric_mean_unweighted() {
        let x = vec![1., 2., 3.];
        let weights = vec![1.0, 1.0, 1.0];
        assert_relative_eq!(weighted_geometric_mean(&x, &weights), 1.8171205928321397);
    }

    #[test]
    fn test_weighted_geometric_mean_weighted() {
        let x = vec![1., 2., 3.];
        let weights = vec![3.0, 2.0, 1.0];
        assert_relative_eq!(weighted_geometric_mean(&x, &weights), 1.5130857494229015);
    }

    #[test]
    fn test_arithmetic_mean() {
        let x = vec![1., 2., 3.];
        assert_relative_eq!(arithmetic_mean(&x), 2.0);
    }
}
