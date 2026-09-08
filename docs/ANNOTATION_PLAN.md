# Siganus fuscescens (Sfus_1.0 / JG2022) — Re-annotation Plan

Goal: produce a community-grade, citable structural + functional annotation of
GCA_042847805.1 (Sfus_1.0), released as GFF3 + GTF + FASTA products with a DOI,
linked from the Kishaz/Siganus_fuscescens-Genome_analysis_pipeline repo.

## 0. Baseline (2024 paper, Mar Biotechnol, doi:10.1007/s10126-024-10325-9)

| Item | Published v1 |
|---|---|
| Assembly | 164 contigs, 498,171,360 bp, N50 7.16 Mb, GC 42.9% |
| Repeats | 17.43% reported |
| Gene models | 28,351 genes, BRAKER3, **protein evidence only** (OrthoDB v11 Vertebrata) |
| RNA-seq evidence | none |
| UTRs | none |
| Functional annotation | none released |
| ncRNA | tRNAscan script present but empty (`14-tRNA-scan.sh` is 0 bytes) |
| Public GFF/GTF | none — assembly on NCBI has **no annotation** |

## 1. Verified facts (checked this session)

- `Sfus-JG2022-Masked.fasta` = 164 seqs, 498,171,360 bp — **byte-exact length match**
  to GCA_042847805.1, contig-for-contig, in the same order. Mapping written to
  `_meta/contig_name_map.tsv` (`ctg000580_np1212` -> `contig_001` -> `BAAFLC010000001.1`).
- Only 5 N bases in the whole assembly (gapless contigs).
- **Softmasked fraction of the supplied FASTA is 5.98%, not the 17.43% in the paper.**
  This file looks like a single masking round, not the cumulative iterative mask.
  Annotating on an under-masked genome is a likely source of the inflated 28,351
  gene count (TE-derived ORFs called as genes). => repeat annotation gets redone.
- RNA-seq available publicly:
  - *S. fuscescens*: **1 run**, SRR36642802 (PRJNA1395358), 7.9 Gbp, paired, NovaSeq X.
  - *S. canaliculatus* (congener, frequently confused/synonymised with S. fuscescens):
    **80 RNA-Seq runs**. Best panel = PRJNA1078125, 8 runs / 164.5 Gbp, multi-tissue
    (SigCg, Ch, Ck, Cl, Cm, Cs, Csp, Ct). Also PRJNA1025886 (M/F gonad), PRJNA774566.
- Congeners with public chromosome-level assemblies for synteny/QC:
  S. canaliculatus (GCA_053572195.1, annotated), S. javus, S. virgatus.

## 2. Coordinate system decision

Primary release uses **GenBank accessions** (`BAAFLC010000001.1` …) so the GTF drops
straight onto the NCBI-downloaded FASTA. Shipped alongside:
- `Sfus_1.0.genomic.fna` (GenBank names, softmasked, produced by us)
- `contig_name_map.tsv` (3-way map)
- convenience GFF3/GTF renamed to `contig_001..164` and to the original `ctg*_np1212`

## 3. Pipeline

### 3.1 Genome preparation
- Rename headers to GenBank accessions; strip description fields.
- `samtools faidx`; record md5 per sequence.

### 3.2 Repeat annotation (redo)
- RepeatModeler2 (`-LTRStruct`) de novo library.
- Merge with Dfam 3.8+ curated Actinopterygii/Teleostei + RepBase-free curated sets.
- Classify unknowns (DeepTE / TEsorter) so the repeat GFF is usable, not "Unknown".
- RepeatMasker `-xsmall` -> cumulative softmasked genome + `repeats.gff3` + `.tbl`.
- Alternative single-command route: **EarlGrey** (wraps the above, better curation).
- Deliverable: honest repeat landscape table, softmasked FASTA at the true %.

### 3.3 Transcript evidence
- Fetch SRR36642802 (S. fuscescens) + selected S. canaliculatus multi-tissue runs.
- `fastp` trim -> `HISAT2` (low-memory; STAR would exceed 30 GB RAM here) to the
  masked genome. Record mapping rate per library; drop libraries < ~70% (cross-species).
- `StringTie` per library + merge -> evidence GTF (also gives UTRs later).

### 3.4 Protein evidence
- OrthoDB v12 Actinopterygii + Vertebrata.
- Curated percomorph proteomes: Danio rerio, Oryzias latipes, Gasterosteus aculeatus,
  Takifugu rubripes, Lates calcarifer, Larimichthys crocea, Sparus aurata,
  Dicentrarchus labrax, **Siganus canaliculatus (GCA_053572195.1)**.
- `miniprot` spliced alignment as an independent evidence track.

### 3.5 Gene prediction
- **BRAKER3** (GeneMark-ETP + AUGUSTUS + TSEBRA) with RNA-seq **and** protein hints.
- **Helixer** (evidence-free deep learning) as an independent second prediction.
- **GALBA/miniprot-based** models for loci the others miss.
- Combine with TSEBRA / AGAT; keep BRAKER3 as the backbone.

### 3.6 Filtering and refinement
- Drop models overlapping repeats above a threshold and lacking protein/RNA support.
- Drop single-exon models with no homology and no expression.
- Add UTRs from StringTie assemblies (GUSHR or direct transfer).
- Enforce: valid start/stop, no internal stops, canonical splice sites where possible.
- Assign stable IDs: `SFUS_G000001` / `SFUS_G000001.1` (documented, versioned).

### 3.7 Non-coding RNA
- tRNAscan-SE 2.0 (finally fills the empty script 14)
- Infernal `cmscan` vs Rfam 15 -> ncRNA GFF3
- barrnap (rRNA), MirMachine (miRNA)

### 3.8 Functional annotation
- InterProScan 5 (Pfam, PANTHER, CDD, SUPERFAMILY, TIGRFAM, SignalP/TMHMM if licensed)
- eggNOG-mapper v2 -> GO, KEGG KO/pathway, COG
- DIAMOND blastp vs UniProt/Swiss-Prot + zebrafish RefSeq -> gene symbols & names
- Merge into GFF3 attributes (`product=`, `Ontology_term=`, `Dbxref=`) via AGAT.
- Revisit the paper's LC-PUFA (elovl/fads) and venom gene stories with real models.

### 3.9 QC / validation
- BUSCO actinopterygii_odb10 **protein mode** on the predicted proteome (target: >95% C)
- OMArk (proteome consistency / contamination / fragmentation)
- `gt gff3validator` + AGAT sanity checks; GTF round-trip test
- Compare gene count, mean exons/gene, mean CDS length, intron length distribution
  against S. canaliculatus, D. labrax, L. crocea. Expect ~22,000–26,000 genes,
  i.e. a *drop* from 28,351 that we can defend.
- Reciprocal-best-hit / synteny check vs annotated S. canaliculatus.

### 3.10 Release
Directory `release/Sfus_1.0.annotation.v2/`:
```
Sfus_1.0.genomic.fna(.gz)          softmasked, GenBank names
Sfus_1.0.annotation.v2.gff3.gz     canonical
Sfus_1.0.annotation.v2.gtf.gz      <- what the requester asked for
Sfus_1.0.proteins.faa.gz
Sfus_1.0.cds.fna.gz
Sfus_1.0.transcripts.fna.gz
Sfus_1.0.ncRNA.gff3.gz
Sfus_1.0.repeats.gff3.gz  + repeat_landscape.tbl
Sfus_1.0.functional_annotation.tsv.gz
contig_name_map.tsv
README.md  MANIFEST.md5  CHANGELOG.md  LICENSE (CC-BY-4.0)
methods.md (exact versions, commands, databases, dates)
```
- Push to GitHub (files >100 MB via Release assets or Git LFS, not the tree).
- Mint a **Zenodo DOI** for the release; add citation block to the repo README.
- Optionally submit the annotation to NCBI (GenBank annotation update for
  GCA_042847805.1) so it is discoverable without GitHub.
- Reply to the person who asked with the DOI + direct links.

## 4. Open decisions
1. Compute venue (this WSL box vs HPC).
2. Any in-house / unpublished S. fuscescens RNA-seq?
3. Whether to include cross-species S. canaliculatus RNA-seq (recommended, flagged).
4. Depth: structural-only fast track vs full functional release.
