# Test data with example for haploinsufficiency (HI)

We create test data for HI using the following data.

## Established HI Genes

We use the following two HI genes on forward/reverse strand for core test data.

- gene: APOB
    - chr2:21224301-21266945
    - reverse strand
- gene: COL3A1
    - chr2:189839099-189877472
    - forward strand
- gene: ACVRL1
    - has pathogenic variants in last coding exon
    - chr12:52301202-52317145
    - forward strand

## Established Non-HI genes

- gene: COL16A1
    - chr1:32117848-32169768
    - reverse strand

## Predicted Non-Established Genes

The following genes are predicted to be HI but are not in ClinGen.
We use them for extended test data.

- gene: MFN2
    - pLI: 0.99
    - NOT in DECIPHER HI
    - chr1:12,040,238-12,073,572
    - forward strand
- gene: REV3L
    - is in DECIPHER HI
    - chr6:111,620,234-111,804,432
    - reverse strand
- gene: RERE
    - is in DECIPHER HI
    - chr1:8,412,464-8,877,699
    - reverse strand

## ClinVar Variants

- we extract all ClinVar variants in all genes from above
- also, we build an "empty" clinvar file with no variants
