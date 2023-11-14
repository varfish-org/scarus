#!/usr/bin/env

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

set -x

cd $SCRIPTPATH

rm -f \
    ClinGen_region_curation_list_GRCh37.tsv \
    ClinGen_region_curation_list_GRCh38.tsv \
    ClinGen_gene_curation_list_GRCh37.tsv

wget -q \
    http://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh37.tsv \
    http://ftp.clinicalgenome.org/ClinGen_region_curation_list_GRCh38.tsv \
    http://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv
