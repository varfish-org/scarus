# Test data with example for HI and TS.

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

# Test data with example for triplosensitivity (TS)

We also incorporate the HI genes as this is a handled in ACMG for CN gain.

We create test data for TS using the following data.

- gene: LMNB1
    - gene: HGNC:6637
    - chr5:126112315-126172712
- gene: PLP1
    - gene: HGNC:9086
    - chrX:103031434-103047548

These are the two only (2023-11-14) established TS genes with level 3.

## Established Non-TS genes

- gene: VCX3A
    - gene: HGNC:18159
    - chrX:6451659-6453159

## Established TS Region

- region: ISCA-37446
    - fully contains COMT, CLDN5 (and others)
    - partially contains PI4KA
    - nearby: USP18 (left)
        - HGNC:12616
    - nearby: SNAP29 (right)
        - HGNC:11133

## Established TS Benign Region

- region: ISCA-46459
    - chr11:89786909-89923863
    - contained: UBTFL1, NAALAD2
    - nearby: TRIM49C (left)
        - HGNC:38877
    - nearby: CHORDC1 (right)
        - HGNC:14525

## Additional Genes

The following lincRNA genes are included for testing protein-coding/non-coding gene distinction.

- gene: MACORIS
  - HGNC:53963
