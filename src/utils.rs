use std::collections::HashMap;

pub fn index_mask(needle: &str, haystack: &[String]) -> Vec<usize> {
    haystack
        .iter()
        .enumerate()
        .filter(|(_, target)| target.contains(needle))
        .map(|(i, _)| i)
        .collect()
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
    group_sizes.into_iter().collect()
}
