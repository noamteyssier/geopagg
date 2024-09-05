use std::collections::HashMap;

use itertools::Itertools;

pub fn index_mask(needle: &str, haystack: &[String]) -> Vec<usize> {
    haystack
        .iter()
        .enumerate()
        .filter(|(_, target)| target.contains(needle))
        .map(|(i, _)| i)
        .collect()
}

/// Transforms a vector of values into z-scores
pub fn zscore_transform(values: &[f64]) -> Vec<f64> {
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let variance = values.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / values.len() as f64;
    let std_dev = variance.sqrt();
    values.iter().map(|x| (x - mean) / std_dev).collect()
}

pub fn select_indices<T: Copy>(indices: &[usize], data: &[T]) -> Vec<T> {
    indices.iter().map(|i| data[*i]).collect()
}

/// Calculates the membership sizes of each gene and the occurence of each membership size
pub fn calculate_group_sizes(genes: &[String], unique_genes: &[&String]) -> Vec<(usize, usize)> {
    let mut group_sizes = HashMap::new();
    for gene in unique_genes {
        let gene_indices = index_mask(gene, genes);
        let membership_size = gene_indices.len();
        if let Some(group_size) = group_sizes.get_mut(&membership_size) {
            *group_size += 1;
        } else {
            group_sizes.insert(membership_size, 1);
        }
    }
    group_sizes.into_iter().sorted().collect()
}

#[cfg(test)]
mod testing {
    use itertools::Itertools;

    #[test]
    fn test_index_mask() {
        let haystack = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        let needle = "a";
        let result = super::index_mask(needle, &haystack);
        assert_eq!(result, vec![0]);
    }

    #[test]
    fn test_select_indices() {
        let data = vec!["a", "b", "c"];
        let indices = vec![0, 2];
        let result = super::select_indices(&indices, &data);
        assert_eq!(result, vec!["a", "c"]);
    }

    #[test]
    fn test_calculate_group_sizes() {
        let genes = vec![
            "a".to_string(),
            "a".to_string(),
            "b".to_string(),
            "c".to_string(),
        ];
        let unique_genes = [String::from("a"), String::from("b"), String::from("c")];
        let unique_gene_references = unique_genes.iter().collect::<Vec<_>>();
        let result = super::calculate_group_sizes(&genes, &unique_gene_references);

        assert!(result.iter().contains(&(1, 2)));
        assert!(result.iter().contains(&(2, 1)));
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_calculate_group_sizes_with_duplicates() {
        let genes = vec![
            "a".to_string(),
            "a".to_string(),
            "b".to_string(),
            "b".to_string(),
            "c".to_string(),
        ];
        let unique_genes = [String::from("a"), String::from("b"), String::from("c")];
        let unique_gene_references = unique_genes.iter().collect::<Vec<_>>();
        let result = super::calculate_group_sizes(&genes, &unique_gene_references);

        assert!(result.iter().contains(&(2, 2)));
        assert!(result.iter().contains(&(1, 1)));
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_zscore_transform() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let z = super::zscore_transform(&x);
        let expected = [
            -f64::sqrt(2.0),
            -0.7071067811865475,
            0.0,
            0.7071067811865475,
            f64::sqrt(2.0),
        ];
        for (a, b) in z.iter().zip(expected.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }
}
