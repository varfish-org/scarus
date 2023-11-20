//! Main entry point for Scarus application.

// #![deny(clippy::pedantic)]
#![allow(clippy::must_use_candidate)]
#![allow(clippy::module_name_repetitions)]
// #![warn(missing_docs)]

use clap::{Parser, Subcommand};

pub mod common;
pub mod seqvars;
pub mod strucvars;

/// CLI parser based on clap.
#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "scarus - ACMG rule evaluation",
    long_about = "This tool provides functionality for automated ACMG rule evaluation"
)]
struct Cli {
    /// Commonly used arguments
    #[command(flatten)]
    common: common::Args,

    /// The sub command to run
    #[command(subcommand)]
    command: Commands,
}

/// Enum supporting the parsing of top-level commands.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Subcommand)]
enum Commands {
    /// Structural variant related commands.
    Strucvars(strucvars::Args),
    /// Sequence variant related commands.
    Seqvars(seqvars::Args),
}

#[tokio::main]
async fn main() -> Result<(), anyhow::Error> {
    let cli = Cli::parse();

    // Build a tracing subscriber according to the configuration in `cli.common`.
    let collector = tracing_subscriber::fmt()
        .with_target(false)
        .with_writer(std::io::stderr)
        .with_max_level(match cli.common.verbose.log_level() {
            Some(level) => match level {
                log::Level::Error => tracing::Level::ERROR,
                log::Level::Warn => tracing::Level::WARN,
                log::Level::Info => tracing::Level::INFO,
                log::Level::Debug => tracing::Level::DEBUG,
                log::Level::Trace => tracing::Level::TRACE,
            },
            None => tracing::Level::INFO,
        })
        .compact()
        .finish();
    tracing::subscriber::set_global_default(collector)?;

    tracing::info!("Starting SCARUS -- taking a chunk out of your variants...");

    match &cli.command {
        Commands::Strucvars(args) => strucvars::run(&cli.common, args)?,
        Commands::Seqvars(args) => seqvars::run(&cli.common, args)?,
    }

    tracing::info!("All done. Have a nice day!");

    Ok(())
}
