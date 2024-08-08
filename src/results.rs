#[derive(Debug, Clone)]
pub struct GeneResult {
    pub gene: String,
    pub wgm: f64,
    pub logfc: f64,
    pub amalgam: bool,
    pub empirical_fdr: f64,
    pub adjusted_empirical_fdr: f64,
}
impl GeneResult {
    /// Creates a new GeneResult
    ///
    /// `empirical_fdr` and `adjusted_empirical_fdr` are initialized to 1.0 and adjusted later in `empirical_fdr` calculations.
    pub fn new(gene: String, wgm: f64, logfc: f64, amalgam: bool) -> Self {
        Self {
            gene,
            wgm,
            logfc,
            amalgam,
            empirical_fdr: 1.0,          // initialized to 1.0
            adjusted_empirical_fdr: 1.0, // initialized to 1.0
        }
    }
}

pub struct GeoPAGGResults {
    pub genes: Vec<String>,
    pub wgms: Vec<f64>,
    pub logfcs: Vec<f64>,
    pub amalgams: Vec<bool>,
    pub empirical_fdr: Vec<f64>,
    pub adjusted_empirical_fdr: Vec<f64>,
}
impl GeoPAGGResults {
    pub fn from_vec(mut gene_results: Vec<GeneResult>) -> Self {
        let mut genes = vec![];
        let mut wgms = vec![];
        let mut logfcs = vec![];
        let mut amalgams = vec![];
        let mut empirical_fdr = vec![];
        let mut adjusted_empirical_fdr = vec![];

        // Sort the results by gene name
        gene_results.sort_unstable_by(|a, b| a.gene.cmp(&b.gene));

        for gene_result in gene_results {
            genes.push(gene_result.gene);
            wgms.push(gene_result.wgm);
            logfcs.push(gene_result.logfc);
            amalgams.push(gene_result.amalgam);
            empirical_fdr.push(gene_result.empirical_fdr);
            adjusted_empirical_fdr.push(gene_result.adjusted_empirical_fdr);
        }

        Self {
            genes,
            wgms,
            logfcs,
            amalgams,
            empirical_fdr,
            adjusted_empirical_fdr,
        }
    }

    pub fn pprint(&self) {
        println!("Gene\tWGM\tLogFC\tAmalgam\tEmpirical FDR\tAdjusted Empirical FDR");
        for i in 0..self.genes.len() {
            println!(
                "{}\t{}\t{}\t{}\t{}\t{}",
                self.genes[i],
                self.wgms[i],
                self.logfcs[i],
                self.amalgams[i],
                self.empirical_fdr[i],
                self.adjusted_empirical_fdr[i]
            );
        }
    }
}
