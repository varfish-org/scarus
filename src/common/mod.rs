//! Commonly used code.

use clap::Parser;
use clap_verbosity_flag::{InfoLevel, Verbosity};

/// Commonly used command line arguments.
#[derive(Parser, Debug)]
pub struct Args {
    /// Verbosity of the program
    #[clap(flatten)]
    pub verbose: Verbosity<InfoLevel>,
}

/// Assembly to be passed on the command line.
#[derive(
    Debug,
    Clone,
    Copy,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    clap::ValueEnum,
    serde::Deserialize,
    serde::Serialize,
)]
pub enum Assembly {
    /// GRCh37
    Grch37,
    /// GRCh38
    Grch38,
}

impl From<Assembly> for biocommons_bioutils::assemblies::Assembly {
    fn from(val: Assembly) -> Self {
        match val {
            // GRCh37p10 is the first one with chrMT.
            Assembly::Grch37 => biocommons_bioutils::assemblies::Assembly::Grch37p10,
            // All canonical contigs including chrMT are in GRCh38.
            Assembly::Grch38 => biocommons_bioutils::assemblies::Assembly::Grch38,
        }
    }
}

impl From<biocommons_bioutils::assemblies::Assembly> for Assembly {
    fn from(val: biocommons_bioutils::assemblies::Assembly) -> Self {
        match val {
            biocommons_bioutils::assemblies::Assembly::Grch37 => Assembly::Grch37,
            biocommons_bioutils::assemblies::Assembly::Grch37p10 => Assembly::Grch37,
            biocommons_bioutils::assemblies::Assembly::Grch38 => Assembly::Grch38,
        }
    }
}
