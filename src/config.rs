use adjustp::{adjust, Procedure};

/// Configurations for weighting p-values in the GeoPAGG algorithm
#[derive(Debug, Clone, Copy)]
pub enum WeightConfig {
    /// Equal weights for all p-values
    Balanced,
    /// Weights increase linearly with rank
    ///
    /// n: number of p-values
    /// i: rank of the p-value (1-indexed)
    ///
    /// w_i = i / n
    RankOrder,
    /// All weights are 1.0 except the first, which is set to alpha
    DropFirst { alpha: f64 },
}

impl WeightConfig {
    /// Builds a vector of weights based on the chosen configuration
    ///
    /// # Arguments
    ///
    /// * `n` - The number of weights to generate
    ///
    /// # Returns
    ///
    /// A vector of f64 weights
    pub fn build_weights(&self, n: usize) -> Vec<f64> {
        match self {
            WeightConfig::Balanced => vec![1.0; n],
            WeightConfig::RankOrder => (1..n + 1).map(|i| (i as f64) / (n as f64)).collect(),
            WeightConfig::DropFirst { alpha } => {
                let mut weights = vec![1.0; n];
                weights[0] = *alpha;
                weights
            }
        }
    }
}

/// Configurations for transforming p-values in the GeoPAGG algorithm
#[derive(Debug, Clone, Copy)]
pub enum TransformConfig {
    /// No transformation
    Identity,
    /// False Discovery Rate correction (Benjamini-Hochberg)
    Fdr,
    /// Bonferroni correction
    Bonferroni,
}

impl TransformConfig {
    /// Transforms a slice of p-values based on the chosen configuration
    ///
    /// # Arguments
    ///
    /// * `pvalues` - A slice of f64 p-values to transform
    ///
    /// # Returns
    ///
    /// A vector of transformed f64 p-values
    #[must_use]
    pub fn transform(&self, pvalues: &[f64]) -> Vec<f64> {
        match self {
            TransformConfig::Identity => pvalues.to_vec(),
            TransformConfig::Fdr => adjust(pvalues, Procedure::BenjaminiHochberg),
            TransformConfig::Bonferroni => adjust(pvalues, Procedure::Bonferroni),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_build_weights_balanced() {
        let weights = WeightConfig::Balanced.build_weights(5);
        assert_eq!(weights, vec![1.0; 5]);
    }

    #[test]
    fn test_build_weights_rank_order() {
        let weights = WeightConfig::RankOrder.build_weights(5);
        assert_eq!(weights, vec![0.2, 0.4, 0.6, 0.8, 1.0]);
    }

    #[test]
    fn test_build_weights_drop_first() {
        let weights = WeightConfig::DropFirst { alpha: 0.5 }.build_weights(5);
        assert_eq!(weights, vec![0.5, 1.0, 1.0, 1.0, 1.0]);
    }
}
