#!/usr/bin/bash

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

set -x
set -euo pipefail

# See README.md for a description of genes and methodology.

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# -- transcripts --

gene_symbols="
ACVRL1
APOB
COL3A1

LMNB1
PLP1

NPHP1
MALL
MTLN
LIMS4

MFN2
REV3L
RERE

MACORIS
"

# wget -O /tmp/cdot-0.2.22.GCF_000001405.25_GRCh37.p13_genomic.105.20201022.gff.json.gz \
#     https://github.com/SACGF/cdot/releases/download/data_v0.2.22/cdot-0.2.22.GCF_000001405.25_GRCh37.p13_genomic.105.20201022.gff.json.gz

rm -rf $SCRIPTPATH/txs_example_hi.bin.zst
mehari db create \
    --genome-release grch37 \
    --path-mane-txs-tsv $SCRIPTPATH/tx-mane.tsv \
    --path-out $SCRIPTPATH/txs_example_hi.bin.zst \
    --path-cdot-json /tmp/cdot-0.2.22.GCF_000001405.25_GRCh37.p13_genomic.105.20201022.gff.json.gz \
    --path-seqrepo-instance ~/mehari-data-tx/seqrepo/master \
    $(for gene_symbol in $gene_symbols; do echo --gene-symbols $gene_symbol; done)

# -- clinvar (filled) --

wget -O $TMPDIR/annonars-clinvar-minimal-grch37-20231112+0.25.0.tar.gz \
  https://github.com/bihealth/annonars-data-clinvar/releases/download/annonars-data-clinvar-20231112/annonars-clinvar-minimal-grch37-20231112+0.25.0.tar.gz
tar -C $TMPDIR -xf $TMPDIR/annonars-clinvar-minimal-grch37-20231112+0.25.0.tar.gz

cat <<EOF | tr ' ' '\t' > $TMPDIR/regions.bed
1 8412464 8877699
1 12040238 12073572
2 21224301 21266945
2 189839099 189877472
6 111620234 111804432
12 52314487 52314714
EOF

mkdir -p tests/data/strucvars/hi_ts/clinvar
rm -rf tests/data/strucvars/hi_ts/clinvar/rocksdb

annonars db-utils copy \
    --path-in $TMPDIR/annonars-clinvar-minimal-grch37-20231112+0.25.0/rocksdb \
    --path-out tests/data/strucvars/hi_ts/clinvar/rocksdb \
    --path-beds $TMPDIR/regions.bed

# -- clinvar (empty) --

mkdir -p tests/data/strucvars/hi_ts/clinvar-empty
rm -rf tests/data/strucvars/hi_ts/clinvar-empty/rocksdb

annonars db-utils copy \
    --path-in $TMPDIR/annonars-clinvar-minimal-grch37-20231112+0.25.0/rocksdb \
    --path-out tests/data/strucvars/hi_ts/clinvar-empty/rocksdb \
    --path-beds /dev/null
