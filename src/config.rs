use adjustp::{adjust, Procedure};

#[derive(Debug, Clone, Copy)]
pub enum WeightConfig {
    Balanced,
    RankOrder,
    DropFirst { alpha: f64 },
}

impl WeightConfig {
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

#[derive(Debug, Clone, Copy)]
pub enum TransformConfig {
    Identity,
    Fdr,
    Bonferroni,
}
impl TransformConfig {
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
