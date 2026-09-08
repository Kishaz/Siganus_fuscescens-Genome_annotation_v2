#!/usr/bin/env bash
# Reachability matrix for every external resource the pipeline needs.
# Run on the login node, then again inside an salloc shell for compute nodes.
echo "== proxy env =="; env | grep -iE 'proxy|no_proxy' || echo "(none set)"
echo "== DNS =="
for h in data.orthodb.org www.dfam.org bioinf.uni-greifswald.de ftp.ebi.ac.uk; do
  printf "%-32s " "$h"; getent hosts "$h" | head -1 || echo "NO-DNS"
done
echo "== HTTP(S) reachability =="
probe(){ printf "%-58s " "$1"; c=$(curl -sS -o /dev/null -w "%{http_code}" --max-time 20 -L "$1" 2>/dev/null); echo "${c:-CONNFAIL}"; }
probe https://ftp.ncbi.nlm.nih.gov/genomes/
probe https://conda.anaconda.org/bioconda/
probe https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/
probe https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/
probe http://eggnog6.embl.de/download/
probe https://busco-data.ezlab.org/v5/data/lineages/
probe https://rest.uniprot.org/
probe https://api.ncbi.nlm.nih.gov/datasets/v2alpha/
probe https://github.com
probe https://objects.githubusercontent.com
probe https://data.orthodb.org/current/download/
probe https://www.dfam.org/releases/current/families/
probe https://dfam.org/releases/
probe https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/
probe https://zenodo.org/api/records
echo "== verbose diagnosis of one failure =="
curl -v --max-time 20 https://data.orthodb.org/ 2>&1 | grep -E '^\*|HTTP/' | head -12
