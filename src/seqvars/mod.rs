//! ACMG rule prediction for sequence variants.

use clap::Parser;

/// Command line arguments for `seqvars` command.
#[derive(Parser, Debug)]
#[command(about = "ACMG rule applications for sequence variants", long_about = None)]
pub struct Args {}

/// Main entry point for the `seqvars` command.
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
    Ok(())
}
