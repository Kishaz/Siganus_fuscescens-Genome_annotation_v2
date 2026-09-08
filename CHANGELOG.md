# Changelog

## v2.0 — pipeline established

First release of the re-annotation pipeline. Supersedes the annotation described
in Mwamburi *et al.* 2024, which was produced by `13-gene prediction.sh` in
[Siganus_fuscescens-Genome_analysis_pipeline](https://github.com/Kishaz/Siganus_fuscescens-Genome_analysis_pipeline)
and whose output files were never released.

### Added

- **Transcript evidence.** BRAKER3 now runs with RNA-seq *and* protein hints;
  v1 was protein-only. Conspecific data only (SRR36642802).
- **Helixer v0.3.7** as a second, evidence-free prediction on GPU, consolidated
  with BRAKER3 by **TSEBRA**.
- **Liftoff** transfer from *Siganus canaliculatus* GCA_053572195.1 — a
  chromosome-level annotated congener that did not exist in 2024 — used as an
  independent evidence axis and a recall check.
- **Evidence-based filtering** on four independent axes (homology, expression,
  congener support, repeat overlap). Models are dropped only when structurally
  invalid or failing every axis. `evidence_support.tsv` ships with the release
  so users can re-filter without rerunning anything.
- **UTRs**, transferred from StringTie assemblies under strict compatibility
  rules (splice-chain agreement, no reach across a neighbouring CDS).
- **Stable identifiers** `SFUS_G000001.1`, position-sorted, with old IDs
  preserved in `old_locus_tag` and a full old→new map.
- **Non-coding RNA**: tRNAscan-SE 2.0.13 with high-confidence filtering,
  barrnap, Infernal/Rfam, and MirMachine. v1 shipped none — its
  `14-tRNA-scan.sh` was a zero-byte file.
- **Functional annotation**: InterProScan, eggNOG-mapper 2.1.15, Swiss-Prot.
  None was released with v1.
- **compleasm** alongside BUSCO 6.1.0, and completeness reported as *recovery of
  the assembly's own BUSCO content* rather than as a bare percentage.
- **OMArk** proteome consistency checking.
- **Acceptance gates** that fail the run rather than shipping a degraded release.
- **Unit tests** for every helper that makes a judgement call.
- `CITATION.cff` and `.zenodo.json` for citable, DOI-minted releases.

### Changed

- **Repeat annotation**: EarlGrey 7.3.1 replaces bare RepeatModeler2 +
  RepeatMasker. v1 left a large `Unknown` fraction, which makes a repeat track
  hard to reuse and makes it impossible to judge which gene models are
  TE-derived. `REPEAT_METHOD=manual` restores the v1 route.
- **Dfam 3.7 → 4.0**, and the new partitioned FamDB format.
- **OrthoDB v11 Vertebrata → v12, Actinopterygii partition.** v12 ships no
  ray-finned partition, so it is cut from Vertebrata by NCBI taxid:
  19,393,872 proteins across 1,006 species → 6,171,398 across 237, a 3.1×
  reduction. Removes tetrapod sequences that are weak splice-site donors for a
  fish while dominating BRAKER's runtime.
- **BRAKER3 → v3.1.1**, pinned by explicit container tag rather than `:latest`.
- **Sequence naming**: releases are on GenBank accessions by default, with
  copies under NCBI sequence names and the original assembly contig names.
- **Portability**: one `config.sh` replaces the `/path/to/your/...` placeholders
  throughout v1. No `#SBATCH --account`/`--partition` is hard-coded.
- **Orchestration**: stages submit as a Slurm DAG, so Helixer runs on GPU while
  BRAKER works on CPU and Liftoff/ncRNA run in parallel.

### Fixed

- The circulated `Sfus-JG2022-Masked.fasta` is only **5.98% softmasked**,
  against **21.15%** for the copy NCBI distributes — it is an under-masked
  intermediate. The pipeline now builds from the GenBank assembly with md5
  verification, and uses 21.15% as an independent QC anchor.
- Infernal `cmscan --fmt 2` column indices: score is field 17 and the overlap
  flag field 20. An earlier draft read fields 15 and 17, silently discarding
  every Rfam hit.
- `busco-data.ezlab.org` serves HTTP 200 with an empty body for directory
  listings; lineage tarballs are now resolved through `file_versions.tsv`,
  which also supplies md5s.
- A `[ -n "$VAR" ] && ...` guard under `set -e` aborted the mirror script
  silently, so Rfam was never fetched despite a clean exit status.

### Known limitations

- Only one conspecific RNA-seq library is public (7.9 Gbp, single tissue).
  Expect roughly 60–70% of loci to carry direct transcript support; the rest
  rest on homology, congener transfer and Helixer. Per-model support is in
  `evidence_support.tsv` — read it before treating a model as confirmed.
- The 80 public *S. canaliculatus* RNA-seq runs are deliberately **not** used as
  transcript evidence, to keep provenance unambiguous. The congener contributes
  only through Liftoff and the protein panel.
- No *Siganus* is present in OrthoDB, so every protein hint is cross-genus.
- Liftoff inherits the reference annotation's errors and cannot see genes the
  reference lacks; it is used as supporting evidence, never to define models.
