#[derive(thiserror::Error, Debug)]
pub enum Error {
    #[error("Unable to find any genes matching substring: {0}")]
    ZeroMatchingGenesFromToken(String),
    #[error("Found a membership size of zero")]
    MembershipSizeOfZero,
    #[error("Found a pvalues len of zero")]
    PvaluesLenZero,
}
