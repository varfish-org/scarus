//! ACMG rule prediction for structural variants.

use clap::Parser;

/// Command line arguments for `strucvars` command.
#[derive(Parser, Debug)]
#[command(about = "ACMG rule applications for structural variants", long_about = None)]
pub struct Args {}

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
    Ok(())
}
