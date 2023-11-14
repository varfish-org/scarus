//! Implementation of evaluation of copy number loss section 1.

use hgvs::data::interface::Provider;
use itertools::Itertools;

use super::result::{Section, L1, L1A, L1B};
use crate::strucvars::{ds::StructuralVariant, eval::del::result::GeneOverlap};

/// Evaluation of deletions, loss of copy number.
///
/// This is mainly used to encapsulate the functionality.  Creating new such
/// objects is very straightforward and cheap.
pub struct Evaluator<'a> {
    /// The parent evaluator.
    parent: &'a super::super::Evaluator,
}

impl<'a> Evaluator<'a> {
    /// Create a new `Evaluator`.
    pub fn with_parent(parent: &'a super::super::Evaluator) -> Self {
        Self { parent }
    }

    /// Perform the evaluation of copy number loss Section 1 and all subsection.
    ///
    /// # Arguments
    ///
    /// * `strucvar` - Structural variant to be evaluated.
    ///
    /// # Returns
    ///
    /// Returns the evaluation results of the section
    ///
    /// # Errors
    ///
    /// If anything goes wrong, it returns a generic `anyhow::Error`.
    pub fn evaluate(&self, strucvar: &StructuralVariant) -> Result<Section, anyhow::Error> {
        let genes = self.overlaps_with_elements(strucvar).map_err(|e| {
            anyhow::anyhow!("issue with overlap computation of {:?}: {}", strucvar, e)
        })?;

        if !genes.is_empty() {
            Ok(Section::L1(L1::L1A(L1A { genes })))
        } else {
            Ok(Section::L1(L1::L1B(L1B::default())))
        }
    }

    /// Determines whether `strucvar` overlaps with any functionally important elements.
    ///
    /// # Arguments
    ///
    /// * `strucvar` - Structural variant to be evaluated.
    ///
    /// # Returns
    ///
    /// Overlapping gene / transcript information.
    ///
    /// # Errors
    ///
    /// When the chromosome name could not be resolved or there was a problem with
    /// accessing the transcript database.
    fn overlaps_with_elements(
        &self,
        strucvar: &StructuralVariant,
    ) -> Result<Vec<GeneOverlap>, anyhow::Error> {
        // Map chromosome name (e.g., chr1) to chromosome accession in this assembly.
        let chrom_acc = self
            .parent
            .chrom_to_ac
            .get(&strucvar.chrom)
            .ok_or_else(|| {
                anyhow::anyhow!("could not resolve chromosome name `{}`", strucvar.chrom)
            })?;

        // Obtain the overlapping transcripts for the given region.
        let txs = self
            .parent
            .provider
            .get_tx_for_region(
                chrom_acc,
                "splign",
                strucvar
                    .start
                    .try_into()
                    .map_err(|e| anyhow::anyhow!("could not convert start position: {}", e))?,
                strucvar
                    .stop
                    .try_into()
                    .map_err(|e| anyhow::anyhow!("could not convert stop position: {}", e))?,
            )
            .map_err(|e| anyhow::anyhow!("problem query transcript database with range: {}", e))?;

        // Extract HGNC ids / transcript ids of coding transcripts.
        let mut gene_txs = txs
            .into_iter()
            .filter_map(|tx| {
                let tx_info = self
                    .parent
                    .provider
                    .get_tx_info(&tx.tx_ac, &tx.alt_ac, "splign")
                    .expect("no tx info?");
                assert!(tx_info.hgnc.starts_with("HGNC:"));
                if tx_info.cds_start_i.is_some() && tx_info.cds_end_i.is_some() {
                    Some((tx_info.hgnc, tx_info.tx_ac))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        gene_txs.sort();

        // Group by and collect for each gene HGNC ID.
        let gene_ovls = gene_txs
            .into_iter()
            .group_by(|(hgnc, _)| hgnc.clone())
            .into_iter()
            .map(|(hgnc, group)| {
                let txs = group.map(|(_, tx)| tx).collect::<Vec<_>>();
                let gene = self
                    .parent
                    .gene_id_data
                    .by_hgnc_id(&hgnc)
                    .expect("could not resolve HGNC ID")
                    .clone();
                GeneOverlap::new(gene, txs)
            })
            .collect::<Vec<_>>();

        tracing::trace!(
            "found gene overlaps for strucvar {:?}: {:?}",
            &strucvar,
            &gene_ovls
        );

        Ok(gene_ovls)
    }
}
