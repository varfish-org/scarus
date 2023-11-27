#!/usr/bin/bash

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

set -x
set -euo pipefail

# See README.md for a description of genes and methodology.

# -- mehari transcripts -----------------------------------------------------

if [[ ! -e $SCRIPTPATH/mehari/txs_example_hi.bin.zst.ok ]]; then
    # The gene symbols to limit the construction to.
    gene_symbols="
    ACVRL1
    APOB
    COL3A1

    LMNB1
    PLP1

    UBTFL1
    NAALAD2
    TRIM49C
    CHORDC1

    MFN2
    REV3L
    RERE

    MACORIS
    "

    # Cleanup any previous output data.
    rm -f $SCRIPTPATH/txs_example_hi.bin.zst $SCRIPTPATH/txs_example_hi.bin.zst.report

    # Download raw data unless a complete download exists.
    if [[ ! -e /tmp/cdot-0.2.22.GCF_000001405.25_GRCh37.p13_genomic.105.20201022.gff.json.gz.ok ]]; then
        wget -O /tmp/cdot-0.2.22.GCF_000001405.25_GRCh37.p13_genomic.105.20201022.gff.json.gz \
            https://github.com/SACGF/cdot/releases/download/data_v0.2.22/cdot-0.2.22.GCF_000001405.25_GRCh37.p13_genomic.105.20201022.gff.json.gz
        touch /tmp/cdot-0.2.22.GCF_000001405.25_GRCh37.p13_genomic.105.20201022.gff.json.gz.ok
    fi

    # Actually run database construction with the given gene symbols.
    mehari db create \
        --genome-release grch37 \
        --path-mane-txs-tsv $SCRIPTPATH/mehari/tx-mane.tsv \
        --path-out $SCRIPTPATH/txs_example_hi.bin.zst \
        --path-cdot-json /tmp/cdot-0.2.22.GCF_000001405.25_GRCh37.p13_genomic.105.20201022.gff.json.gz \
        --path-seqrepo-instance ~/mehari-data-tx/seqrepo/main \
        $(for gene_symbol in $gene_symbols; do echo --gene-symbols $gene_symbol; done)

    # Touch the output marker file.
    touch $SCRIPTPATH/mehari/txs_example_hi.bin.zst.ok
fi

# -- ClinVar seqvars --------------------------------------------------------

if [[ ! -e $SCRIPTPATH/clinvar/rocksdb.ok ]]; then
    # Download raw data unless a complete download exists.
    if [[ ! -e /tmp/annonars-clinvar-minimal-grch37-20231112+0.29.0.ok ]]; then
        wget -O /tmp/annonars-clinvar-minimal-grch37-20231112+0.29.0.tar.gz \
            https://github.com/bihealth/annonars-data-clinvar/releases/download/annonars-data-clinvar-20231112/annonars-clinvar-minimal-grch37-20231112+0.29.0.tar.gz
        tar -C /tmp -xf /tmp/annonars-clinvar-minimal-grch37-20231112+0.29.0.tar.gz
        touch /tmp/annonars-clinvar-minimal-grch37-20231112+0.29.0.ok
    fi

    # Cleanup any previous output data.
    rm -rf $SCRIPTPATH/clinvar/rocksdb

    # Copy selected data.
    annonars db-utils copy \
        --path-in /tmp/annonars-clinvar-minimal-grch37-20231112+0.29.0/rocksdb \
        --path-out $SCRIPTPATH/clinvar/rocksdb \
        --path-beds $SCRIPTPATH/clinvar/regions.bed

    # Touch the output marker file.
    touch $SCRIPTPATH/clinvar/rocksdb.ok
fi

# -- ClinVar seqvars (empty) -------------------------------------------------
#
# This is useful when we want to guarantee that there is no ClinVar data

if [[ ! -e $SCRIPTPATH/clinvar-empty/rocksdb.ok ]]; then
    # Cleanup any previous output data.
    rm -rf $SCRIPTPATH/clinvar-empty/rocksdb

    # Copy null data
    annonars db-utils copy \
        --path-in $SCRIPTPATH/clinvar/rocksdb \
        --path-out $SCRIPTPATH/clinvar-empty/rocksdb \
        --path-beds /dev/null

    # Touch the output marker file.
    touch $SCRIPTPATH/clinvar-empty/rocksdb.ok
fi

# -- ClinVar SV Data --------------------------------------------------------
#
# We import all data as there is not too much of it.

if [[ ! -e $SCRIPTPATH/clinvar-sv/rocksdb.ok ]]; then
    # Download raw data unless a complete download exists.
    if [[ ! -e /tmp/clinvar-data-extract-vars-20231112+0.12.0.ok ]]; then
        wget -O /tmp/clinvar-data-extract-vars-20231112+0.12.0.tar.gz \
            https://github.com/bihealth/clinvar-data-jsonl/releases/download/clinvar-weekly-20231112/clinvar-data-extract-vars-20231112+0.12.0.tar.gz
        tar -C /tmp -xvf /tmp/clinvar-data-extract-vars-20231112+0.12.0.tar.gz
        touch /tmp/clinvar-data-extract-vars-20231112+0.12.0.ok
    fi

    # Cleanup any previous output data
    rm -rf $SCRIPTPATH/clinvar-sv/rocksdb

    # Actually run database construction.
    annonars clinvar-sv import \
        --genome-release grch37 \
        --path-out-rocksdb $SCRIPTPATH/clinvar-sv/rocksdb \
        --path-in-jsonl /tmp/clinvar-data-extract-vars-20231112+0.12.0/clinvar-variants-grch37-strucvars.jsonl.gz

    # Touch the output marker file.
    touch $SCRIPTPATH/clinvar-sv/rocksdb.ok
fi

# -- RefSeq Functional Elements ---------------------------------------------

if [[ ! -e $SCRIPTPATH/functional/rocksdb.ok ]]; then
    # Download raw data unless a complete download exists.
    if [[ ! -e /tmp/GCF_000001405.25_GRCh37.p13_genomic.gff.gz.ok ]]; then
        wget \
            -O /tmp/GCF_000001405.25_GRCh37.p13_genomic.gff.gz \
            https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/9606/105.20201022/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz
        zgrep '^#\|RefSeqFE' \
            /tmp/GCF_000001405.25_GRCh37.p13_genomic.gff.gz \
            > /tmp/GCF_000001405.25_GRCh37.p13_genomic.functional.gff
        touch /tmp/GCF_000001405.25_GRCh37.p13_genomic.functional.gff.ok
    fi

    # Cleanup any previous output data.
    rm -rf $SCRIPTPATH/functional/rocksdb

    # Actually run database construction.
    annonars functional import \
        --genome-release grch37 \
        --path-in-gff /tmp/GCF_000001405.25_GRCh37.p13_genomic.functional.gff \
        --path-out-rocksdb $SCRIPTPATH/functional/rocksdb

    # Touch the output marker file.
    touch $SCRIPTPATH/functional/rocksdb.ok
fi

# -- annonars genes ---------------------------------------------------------

if [[ ! -e $SCRIPTPATH/genes/rocksdb.ok ]]; then
    # Cleanup any previous output data.
    rm -rf $SCRIPTPATH/genes/rocksdb

    # Download the data from S3
    mkdir -p $SCRIPTPATH/genes/rocksdb
    s5cmd \
        --endpoint-url=https://ceph-s3-public.cubi.bihealth.org \
        --no-sign-request \
        sync \
        "s3://varfish-public/full/annonars/genes-3.1+2.1.1+4.4+20230606+10.1+20231123+0.29.3/rocksdb/*" \
        $SCRIPTPATH/genes/rocksdb/

    # Touch the output marker file.
    touch $SCRIPTPATH/genes/rocksdb.ok
fi

# -- annonars regions -------------------------------------------------------

if [[ ! -e $SCRIPTPATH/regions/rocksdb.ok ]]; then
    # Cleanup any previous output data.
    rm -rf $SCRIPTPATH/regions/rocksdb

    # Download the data from S3
    mkdir -p $SCRIPTPATH/regions/rocksdb
    s5cmd \
        --endpoint-url=https://ceph-s3-public.cubi.bihealth.org \
        --no-sign-request \
        sync \
        "s3://varfish-public/full/annonars/regions-grch37-20231122+0.29.3/rocksdb/*" \
        $SCRIPTPATH/regions/rocksdb/

    # Touch the output marker file.
    touch $SCRIPTPATH/regions/rocksdb.ok
fi
