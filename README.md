# *Siganus fuscescens* genome annotation v2

Reproducible structural and functional re-annotation of
[GCA_042847805.1 (Sfus_1.0)](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_042847805.1/),
the *Siganus fuscescens* (Houttuyn, 1782) reference assembly — the first genome
published for the family Siganidae
([Mwamburi *et al.* 2024, *Marine Biotechnology*](https://doi.org/10.1007/s10126-024-10325-9)).

**The assembly is public but carries no annotation record on NCBI.** Anyone who
downloads `GCA_042847805.1` gets sequence and nothing else. This pipeline produces
the GFF3, GTF and derived products the assembly has been missing, in a form that
can be cited and reused without contacting the authors.

Companion to
[Kishaz/Siganus_fuscescens-Genome_analysis_pipeline](https://github.com/Kishaz/Siganus_fuscescens-Genome_analysis_pipeline),
which covers the assembly itself.

---

## What v2 changes

The 2024 annotation ran BRAKER3 with protein evidence only, against OrthoDB v11
Vertebrata, on a genome masked with Dfam 3.7 — and released no annotation files.
v2 rebuilds all of it.

| | v1 (2024) | v2 |
|---|---|---|
| Repeat annotation | RepeatModeler2 + RepeatMasker, Dfam 3.7, large `Unknown` fraction | **EarlGrey 7.3.1** (BLAST-extend-extend curation + classification), **Dfam 4.0** |
| Repeat products | mask only | mask **+ classified `repeats.gff3` + divergence landscape** |
| Transcript evidence | **none** | conspecific RNA-seq → HISAT2 + StringTie 3 |
| Protein evidence | OrthoDB v11 Vertebrata (all 19.4 M) | **OrthoDB v12, Actinopterygii partition** (6.17 M, cut by taxid) + 9 curated percomorph proteomes |
| Gene prediction | BRAKER3 (protein-only) | **BRAKER3 v3.1.1** (RNA-seq **+** protein) |
| Second, independent prediction | — | **Helixer v0.3.7** (deep learning, GPU, evidence-free) |
| Congener evidence | — | **Liftoff** from *S. canaliculatus* GCA_053572195.1 |
| Model consolidation | — | **TSEBRA** |
| Filtering | none — 28,351 genes as emitted | four independent evidence axes; every model ships its support |
| UTRs | none | transferred from assembled transcripts, conservatively |
| Identifiers | `g1234.t1` (run-specific) | `SFUS_G000001.1` (stable, versioned, mapped to old IDs) |
| Non-coding RNA | none (`14-tRNA-scan.sh` was an empty file) | tRNAscan-SE 2.0.13, barrnap, Infernal/Rfam, **MirMachine** |
| Functional annotation | **none released** | InterProScan + eggNOG-mapper 2.1.15 + Swiss-Prot |
| Completeness QC | BUSCO (genome mode only) | **BUSCO 6.1.0 + compleasm**, two lineages, reported as *recovery of the assembly's own BUSCO content* |
| Proteome QC | — | **OMArk** |
| Reproducibility | shell scripts with `/path/to/your/...` placeholders | one `config.sh`, pinned containers, frozen env locks, unit tests |

The gene count is expected to **fall** from 28,351 toward the 20,000–27,000
teleost range. That is the intended outcome, not a regression — v1 had no
transcript evidence and no TE-derived-model filtering. Every dropped model is
retained in `dropped.gff3` with a stated reason.

### Why two gene predictors

BRAKER3 is driven entirely by extrinsic evidence. Helixer uses none — it predicts
from sequence alone. Where they independently agree, the agreement is meaningful;
where BRAKER has no evidence to work with, Helixer still produces a call. With
only one public conspecific RNA-seq library, that gap is real. TSEBRA merges them
with BRAKER3 as the backbone rather than naively unioning them.

### Why a congener axis in filtering

With a single RNA-seq library, homology and expression tend to fail *together*
for one underlying reason: a poorly conserved, tissue-restricted gene fails both.
Counting them as two independent votes overstates the evidence. The Liftoff
transfer from a chromosome-level congener gives such loci a genuinely separate
way to prove they are real. `tests/test_filter.sh` contains the differential case.

---

## Provenance

The pipeline downloads the assembly from GenBank rather than trusting a local
copy, so the release is provably built on the bytes the public gets.

Verified during development:

- The locally held `Sfus-JG2022-Masked.fasta` and `GCA_042847805.1` are
  **base-for-base identical across all 164 sequences** (case-insensitive
  per-sequence md5).
- The contig map (`ctg000580_np1212` → `contig_001` → `BAAFLC010000001.1`) is
  derived from the NCBI assembly report and keyed on unique sequence length, so
  it cannot silently mis-pair.
- The locally circulated FASTA is **5.98% softmasked**; the copy NCBI
  distributes is **21.15%**. The former is an under-masked intermediate. 21.15%
  is used as an independent QC anchor for stage 02.
- **OrthoDB contains no *Siganus*** — every protein hint is cross-genus, which
  independently corroborates this being the first Siganidae genome.

## Requirements

- Slurm cluster, Linux x86-64
- `apptainer` or `singularity` (BRAKER3 and Helixer both run from containers)
- One GPU for stage 07 (optional — set `USE_HELIXER=0` to skip)
- Outbound HTTPS from the login node — see *Firewalled clusters* below if partial
- ~350 GB of working space: ~150 GB databases, ~150 GB intermediates
- No root: everything installs into the project directory via micromamba

Storage is split into two tiers. `SCRATCH` holds regenerable work and databases
and is auto-detected — a candidate must be writable, have ≥350 GB free, and not
be tmpfs, since siting intermediates on a RAM-backed filesystem exhausts node
memory. `PROJ` holds the release, logs and checkpoints. Because scratch is
typically purgeable, expensive artefacts (the softmasked genome, curated repeat
library, and the raw predictions) are checkpointed to `PROJ` as they are
produced, and every consuming stage restores from there if scratch has been
cleared — so a purge costs a copy rather than a rerun.

## Quick start

```bash
git clone <this repo>
cd Siganus_fuscescens-Genome_annotation_v2
$EDITOR config.sh          # SLURM_ACCOUNT, PART_CPU, PART_GPU, PROJ

bash tests/test_filter.sh && bash tests/test_utr_ids.sh && bash tests/test_release.sh

./scripts/00_setup.sh              # envs + containers (~1-2 h)
./scripts/01_prepare_genome.sh     # fetch and verify the assembly
./scripts/03_fetch_evidence.sh     # databases + RNA-seq (login node)

./scripts/preflight.sh             # verify everything BEFORE queuing days of work
DRYRUN=1 ./scripts/run_all.sh      # inspect the submission plan
./scripts/run_all.sh               # submit the DAG
```

### Firewalled clusters

Many HPC systems allowlist only a few external hosts. On the cluster this was
developed against, the reachable set is NCBI, Anaconda, Docker Hub and eggNOG —
while `dfam.org`, `data.orthodb.org`, `bioinf.uni-greifswald.de`,
`busco-data.ezlab.org` and **`ftp.ebi.ac.uk`** all fail to connect, from login
*and* compute nodes alike. EBI being blocked matters most, because InterProScan,
Swiss-Prot and Rfam live there and have no alternative source.

Mirror them from any machine with a working connection (~13.5 GB total):

```bash
./scripts/mirror_fetch.sh /tmp/staging
rsync -avP /tmp/staging/ user@cluster:/path/to/aigo_genome/db/staged/
```

`03_fetch_evidence.sh` prefers `db/staged/` over downloading, and `preflight.sh`
reports exactly which databases are present and what each missing one costs.

No compute stage performs network I/O; all fetching happens on the login node.

## Stages

| # | Stage | Where | Time |
|---|---|---|---|
| 00 | environments, containers, version freeze | login | 1–2 h |
| 01 | fetch + md5-verify assembly | login | 10 min |
| 02 | **EarlGrey** repeats + softmask | cpu | **3–7 d** |
| 03 | RNA-seq, proteomes, databases | login | 2–6 h |
| 04 | fastp → HISAT2 → StringTie 3 | cpu | 4–8 h |
| 05 | OrthoDB clade cut + miniprot | cpu | 2–4 h |
| 06 | **BRAKER3 v3.1.1** (RNA-seq + protein) | cpu | 1–2 d |
| 07 | **Helixer v0.3.7** | **gpu** | 4–8 h |
| 08 | **Liftoff** from *S. canaliculatus* | cpu | 2–4 h |
| 09 | **TSEBRA** consolidation | cpu | 1–2 h |
| 10 | evidence filtering, UTRs, stable IDs | cpu | 1–2 h |
| 11 | tRNAscan-SE / barrnap / Rfam / MirMachine | cpu | 4–8 h |
| 12 | InterProScan + eggNOG + Swiss-Prot | cpu | 1–2 d |
| 13 | BUSCO + compleasm + OMArk + gates | cpu | 2–4 h |
| 14 | package release | login | 15 min |

Every stage is idempotent — rerunning skips completed work.

## Acceptance gates

Stage 13 **fails the run** rather than shipping if:

- BUSCO (`actinopterygii_odb10`, protein mode) completeness < 92%
- gene count outside 20,000–27,000
- mean exons per transcript < 6.0
- any internal stop codons
- mono-exonic transcript fraction > 0.35
- GFF3 does not validate under `gt gff3validator`

`actinopterygii_odb10` is the gated lineage because that is what the 2024 paper
measured — comparing across lineages would make the before/after meaningless.
`odb12` is run and reported alongside.

## Release products

```
Sfus_1.0.annotation.v2.0.gff3.gz     canonical
Sfus_1.0.annotation.v2.0.gtf.gz      the format most people ask for
Sfus_1.0.proteins.faa.gz  .cds.fna.gz  .transcripts.fna.gz
Sfus_1.0.ncRNA.gff3.gz    Sfus_1.0.repeats.gff3.gz
Sfus_1.0.functional_annotation.tsv.gz
Sfus_1.0.evidence_support.tsv.gz     per-model evidence, so users can re-filter
contig_name_map.tsv    alt_names/     GenBank / NCBI / original contig naming
METHODS.md  software_versions.tsv  MANIFEST.md5  LICENSE
```

Coordinates are on **GenBank accessions** by default, so the GTF drops straight
onto the FASTA people download from NCBI. Copies under the other two naming
schemes ship alongside.

## Tests

```bash
bash tests/test_filter.sh     # evidence scoring, incl. the congener differential
bash tests/test_utr_ids.sh    # UTR transfer (both strands), stable IDs
bash tests/test_release.sh    # renaming, ncRNA precedence, functional merge, gates
```

Seconds to run, on synthetic fixtures — no cluster or databases needed.

## Citation

Please cite the genome paper:

> Mwamburi, S.M., Kawato, S., Furukawa, M., Konishi, K., Nozaki, R., Hirono, I.,
> Kondo, H. (2024). De Novo Assembly and Annotation of the *Siganus fuscescens*
> (Houttuyn, 1782) Genome: Marking a Pioneering Advance for the Siganidae Family.
> *Marine Biotechnology*. https://doi.org/10.1007/s10126-024-10325-9

## Licence

Pipeline code MIT; annotation data products CC-BY 4.0. See `LICENSE`.
