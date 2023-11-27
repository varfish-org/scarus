//! ACMG rule prediction for structural variants.

pub mod ds;
pub mod eval;

use std::path::PathBuf;

use clap::Parser;

use crate::{common::Assembly, strucvars::eval::Paths};

use self::ds::StructuralVariant;

/// Command line arguments for `strucvars` command.
#[derive(Parser, Debug)]
#[command(about = "ACMG rule applications for structural variants", long_about = None)]
pub struct Args {
    /// The used assembly.
    #[clap(long, value_enum)]
    pub assembly: Assembly,

    /// Path to Mehari transcript database file.
    #[clap(long)]
    pub path_tx_db: PathBuf,
    /// Path to annonars genes RocksDB.
    #[clap(long)]
    pub path_annonars_genes: PathBuf,
    /// Path to annonars functional elements RocksDB.
    #[clap(long)]
    pub path_annonars_functional: PathBuf,
    /// Path to annonars regions RocksDB.
    #[clap(long)]
    pub path_annonars_regions: PathBuf,

    /// Path to ClinVar seqvars RocksDB.
    #[clap(long)]
    pub path_clinvar_seqvars: PathBuf,
    /// Path to ClinVar strucvars RocksDB.
    #[clap(long)]
    pub path_clinvar_strucvars: PathBuf,
    /// Path to gnomAD SV exomes.
    #[clap(long)]
    pub path_gnomad_sv_exomes: PathBuf,
    /// Path to gnomAD SV genomes.
    #[clap(long)]
    pub path_gnomad_sv_genomes: PathBuf,

    /// The variants to prioritize.
    #[clap(required = true)]
    pub variants: Vec<StructuralVariant>,
}

/// Main entry point for the `strucvars` command.
///
/// # Arguments
///
/// * `common_args` - Commonly used command line arguments.
/// * `args` - Command line arguments specific to `strucvars` command.
///
/// # Errors
///
/// If anything goes wrong, it returns a generic `anyhow::Error`.
pub fn run(common_args: &crate::common::Args, args: &Args) -> Result<(), anyhow::Error> {
    tracing::info!("  running command `strucvars`");
    tracing::info!("  common_args = {:?}", &common_args);
    tracing::info!("  args = {:?}", &args);

    let paths = Paths {
        path_tx_db: args.path_tx_db.clone(),
        path_annonars_genes: args.path_annonars_genes.clone(),
        path_annonars_functional: args.path_annonars_functional.clone(),
        path_annonars_regions: args.path_annonars_regions.clone(),
        path_clinvar_seqvars: args.path_clinvar_seqvars.clone(),
        path_clinvar_strucvars: args.path_clinvar_strucvars.clone(),
        path_gnomad_sv_exomes: args.path_gnomad_sv_exomes.clone(),
        path_gnomad_sv_genomes: args.path_gnomad_sv_genomes.clone(),
    };
    let config = Default::default();

    let evaluator = eval::Evaluator::new(args.assembly.into(), paths, config)
        .map_err(|e| anyhow::anyhow!("failed to create evaluator: {}", e))?;
    for variant in &args.variants {
        tracing::info!("- assessing {:?}", variant);
        let result = evaluator.evaluate(variant)?;
        println!("{}", serde_json::to_string(&result)?);
    }

    Ok(())
}

#[cfg(test)]
mod test {
    #[test]
    fn run_smoke() -> Result<(), anyhow::Error> {
        let common = crate::common::Args {
            verbose: clap_verbosity_flag::Verbosity::new(1, 0),
        };

        let args = super::Args {
            assembly: crate::common::Assembly::Grch37,
            path_tx_db: "tests/data/strucvars/txs_example_hi.bin.zst".into(),
            path_annonars_genes: "tests/data/strucvars/genes/rocksdb".into(),
            path_annonars_functional: "tests/data/strucvars/functional/rocksdb".into(),
            path_annonars_regions: "tests/data/strucvars/regions/rocksdb".into(),
            path_clinvar_seqvars: "tests/data/strucvars/clinvar/rocksdb".into(),
            path_clinvar_strucvars: "tests/data/strucvars/clinvar-sv/rocksdb".into(),
            path_gnomad_sv_exomes: "tests/data/strucvars/gnomad-sv/exac-cnv/rocksdb".into(),
            path_gnomad_sv_genomes: "tests/data/strucvars/gnomad-sv/gnomad-sv2/rocksdb".into(),
            variants: vec!["DEL:1:21000:26000".parse()?, "DUP:1:21000:26000".parse()?],
        };

        super::run(&common, &args)
    }
}
