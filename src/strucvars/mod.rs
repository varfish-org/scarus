//! ACMG rule prediction for structural variants.

pub mod data;
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
    /// Path to HGNC xlink file.
    #[clap(long)]
    pub path_hgnc: PathBuf,
    /// Path to ClinGen dosage genes TSV file.
    #[clap(long)]
    pub path_clingen_dosage_genes: PathBuf,
    /// Path to ClinGen dosage regions TSV file.
    #[clap(long)]
    pub path_clingen_dosage_regions: PathBuf,
    /// Path to Decipher HI predictions TSV file.
    #[clap(long)]
    pub path_decipher_hi: PathBuf,
    /// Path to gnomAD constraints TSV file.
    #[clap(long)]
    pub path_gnomad_constraints: PathBuf,
    /// Path to ClinVar "minimal" RocksDB.
    #[clap(long)]
    pub path_clinvar_minimal: PathBuf,
    /// Path to functional element RocksDB.
    #[clap(long)]
    pub path_functional: PathBuf,
    /// Path to ClinVar SV RocksDB.
    #[clap(long)]
    pub path_clinvar_sv: PathBuf,
    /// Path to gnomAD SV genomes.
    #[clap(long)]
    pub path_gnomad_sv_genomes: PathBuf,
    /// Path to gnomAD SV exomes.
    #[clap(long)]
    pub path_gnomad_sv_exomes: PathBuf,

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
        path_hgnc: args.path_hgnc.clone(),
        path_clingen_dosage_genes: args.path_clingen_dosage_genes.clone(),
        path_clingen_dosage_regions: args.path_clingen_dosage_regions.clone(),
        path_decipher_hi: args.path_decipher_hi.clone(),
        path_gnomad_constraints: args.path_gnomad_constraints.clone(),
        path_clinvar_minimal: args.path_clinvar_minimal.clone(),
        path_functional: args.path_functional.clone(),
        path_clinvar_sv: args.path_clinvar_sv.clone(),
        path_gnomad_sv_genomes: args.path_gnomad_sv_genomes.clone(),
        path_gnomad_sv_exomes: args.path_gnomad_sv_exomes.clone(),
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
            path_tx_db: "tests/data/strucvars/hi_ts/txs_example_hi.bin.zst".into(),
            path_hgnc: "tests/data/hgnc.tsv".into(),
            path_clingen_dosage_regions:
                "tests/data/strucvars/ClinGen_region_curation_list_GRCh37.tsv".into(),
            path_clingen_dosage_genes: "tests/data/strucvars/ClinGen_gene_curation_list_GRCh37.tsv"
                .into(),
            path_decipher_hi: "tests/data/strucvars/decipher_hi_prediction.tsv".into(),
            path_gnomad_constraints: "tests/data/strucvars/gnomad_constraints.tsv".into(),
            path_clinvar_minimal: "tests/data/strucvars/hi_ts/clinvar/rocksdb".into(),
            path_functional: "tests/data/strucvars/hi_ts/functional/rocksdb".into(),
            path_clinvar_sv: "tests/data/strucvars/hi_ts/clinvar-sv/rocksdb".into(),
            path_gnomad_sv_genomes: "tests/data/strucvars/hi_ts/gnomad-sv/gnomad-sv2/rocksdb"
                .into(),
            path_gnomad_sv_exomes: "tests/data/strucvars/hi_ts/gnomad-sv/exac-cnv/rocksdb".into(),
            variants: vec!["DEL:1:21000:26000".parse()?, "DUP:1:21000:26000".parse()?],
        };

        super::run(&common, &args)
    }
}
