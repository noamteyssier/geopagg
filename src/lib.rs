//! GeoPAGG: Geometric P-Value Aggregation for Gene Groups
//!
//! This library implements the GeoPAGG algorithm, which is used for analyzing
//! sgRNA differential abundance data by aggregating p-values and log fold changes
//! across multiple sgRNAs targeting the same gene.
//!
//! The main components of this library are:
//! - `GeoPAGG`: The main algorithm implementation
//! - `WeightConfig`: Configuration for p-value weighting schemes
//! - `TransformConfig`: Configuration for p-value transformation methods
//! - `GeoPAGGResults`: Structure to hold and display the results

mod amalgam;
mod config;
mod geopagg;
mod math;
mod results;
mod utils;

pub use config::{TransformConfig, WeightConfig};
pub use geopagg::GeoPAGG;
pub use results::GeoPAGGResults;
