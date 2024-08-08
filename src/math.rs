use crate::{config::WeightConfig, results::GeneResult};

/// Calculates the weighted geometric mean of a set of values.
///
/// # Arguments
///
/// * `x` - A slice of f64 values to calculate the mean from.
/// * `weights` - A slice of f64 weights corresponding to each value in `x`.
///
/// # Returns
///
/// The calculated weighted geometric mean as an f64.
///
/// # Description
///
/// The weighted geometric mean is calculated as:
/// WGM = exp((Σ(w_i * ln(x_i))) / Σ(w_i))
/// where w_i are the weights and x_i are the values.
///
/// # Panics
///
/// This function will panic if:
/// - The input slices have different lengths.
/// - Any value in `x` is less than or equal to 0.
/// - The sum of weights is 0.
pub fn weighted_geometric_mean(x: &[f64], weights: &[f64]) -> f64 {
    let sum = x
        .iter()
        .zip(weights.iter())
        .map(|(x, w)| x.ln() * w)
        .sum::<f64>();
    let sum_of_weights = weights.iter().sum::<f64>();
    (sum / sum_of_weights).exp()
}

/// Calculates the arithmetic mean of a set of values.
///
/// # Arguments
///
/// * `x` - A slice of f64 values to calculate the mean from.
///
/// # Returns
///
/// The calculated arithmetic mean as an f64.
///
/// # Description
///
/// The arithmetic mean is calculated as:
/// AM = (Σ x_i) / n
/// where x_i are the values and n is the number of values.
///
/// # Panics
///
/// This function will panic if the input slice is empty.
pub fn arithmetic_mean(x: &[f64]) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}

/// Aggregates p-values using a weighted geometric mean approach.
///
/// # Arguments
///
/// * `pvalues` - A mutable slice of f64 p-values to aggregate.
/// * `weight_config` - The configuration for determining weights.
///
/// # Returns
///
/// The aggregated p-value as an f64.
///
/// # Description
///
/// This function performs the following steps:
/// 1. Sorts the p-values in ascending order.
/// 2. Builds weights based on the provided `WeightConfig`.
/// 3. Calculates the weighted geometric mean of the sorted p-values.
///
/// # Panics
///
/// This function will panic if:
/// - Any p-value is less than or equal to 0.
/// - The input slice is empty.
pub fn aggregate_pvalues(pvalues: &mut [f64], weight_config: WeightConfig) -> f64 {
    // Sort the pvalues
    pvalues.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

    let weights = weight_config.build_weights(pvalues.len());
    weighted_geometric_mean(pvalues, &weights)
}

/// Calculates the empirical False Discovery Rate (FDR) for each gene result.
///
/// # Arguments
///
/// * `gene_results` - A mutable slice of `GeneResult` structures.
///
/// # Description
///
/// This function performs the following steps:
/// 1. Sorts the gene results by their weighted geometric mean (WGM) in ascending order.
/// 2. For each gene result:
///    - Counts the number of amalgams with a lower or equal WGM.
///    - Calculates the empirical FDR as (number of amalgams) / (current rank).
///    - Calculates the adjusted empirical FDR as `max(empirical FDR, WGM)`.
///
/// # Note
///
/// This function modifies the input `gene_results` in place, updating the `empirical_fdr`
/// and `adjusted_empirical_fdr` fields of each `GeneResult`.
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
