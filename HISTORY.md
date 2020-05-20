1.0.0 (TBR):

- Integrated [oncoviruses](https://github.com/umccr/oncoviruses). It detects oncoviral content and possible integration sites, as well as genes affected by integration. Reported viral strains in MultiQC and integration sites and affected genes in the cancer report.
- More DRAGEN input options:
    - Support tumor and normal sample names different from the output prefix (read from the --RGSM tags)
    - Support multi-fastq runs (when input fastq files and RGSM are specified in CSV files and put in the DRAGEN outout folder by the UMCCR worklow).
- Custom input. You can also feed custom files using multiple positional arguments. VCF and BAM files are supported. Sample name will be extracted from VCF and BAM headers. For now, VCF is assumed to contain T/N somatic small variant calls, BAM file is assumed to be from tumor. Example:

```
umccrise umccrise sample1.bam sample2.bam sample1.vcf.gz sample3.vcf.gz -o umccrised -j10
```

- Improved parallelization.
    Re-engineered how the number of CPUs is determined, so when we have many cores available, they are not wasted on a single job that doesn't benefit from many threads. Say, PURPLE COBALT optional number is around 10-15, so if we have 28 cores, we will run COBALT and AMBER in parallel with 14 cores each. This should work for local/AWS runs as well as the HPC runs.
- PURPLE updated:
    - Use PURPLE v2.40 that reports tumor mutational load and burden
    - Update COBALT GC profile reference file so it fixes a bug ("We have found an issue with the mappability file which is used by COBALT to decide which windows to fit. The file was created on the full set ref genome, when we should have used the analysis set. The impact is that some large (and important) sections of hg38 which have alt contigs are assigned 0 mappability and filtered by COBALT. Particularly the HLA and ERBB2 regions are affected. We have generated new files already and I have checked them, but we are just testing on some real samples.")
- Integrated/oncoviruses/oncoviral_breakpoints.tsv OptiType for HLA typing (disabled by default, can be enabled with `-T immuno`)
- Re-engineered usage of stages. Add --include-stage (-T) and --exclude-stages (-E), allowed to be provided multiple times. Stages as positional arguments are no longer supported. Instead of:

```
umccrise /bcbio/final/ multiqc structural
```

Use:

```
umccrise /bcbio/final/ -T multiqc -T structural
```

Positional arguments can now be used for extra custom inputs - BAMs or VCFs, e.g.:

```
umccrise sample1.bam sample2.bam sample1.vcf.gz sample3.vcf.gz
```

Sample name will be extracted from VCF and BAM headers. For now, VCF is assumed to contain T/N somatic small variant calls, BAM file is assumed to be from tumor.

In concordance with stages, sample and batch inclusion/exclusion parameters (`-s` and `-e`) are now implemented in the same way allowing to be provided multiple times instead of a comma-separated string. Example:

```
# of all samples in a project, takes only sample1 and sample3, plus all corresponding normal/tumor matches:
umccrise /input/project/final -s sample1 -s sample3

# takes all samples in a project, excluding sample1 and sample2 and corresponding normal/tumor matches:
umccrise /input/project -e sample1 -e sample2
```

Also add stages: `cpsr` (decoupled from `pcgr`), `somatic`, `maf` and `germline` (fine-grained `small_variants`), `mosdepth`, `goleft` and `cacao` (fine-grained `coverage`), `oncoviruses`, `microbiome`, `immuno` (novel stages).

Not that some of the stages are inter-dependent, e.g. `germline` uses `somatic` to get extra germline leakage, `coverage` uses `purple` to get purity-adjusted coverage thresholds, `multiqc` reports results of  `conpair`, `somatic`, `germline` and `oncoviral` stages, `pcgr` uses `purple` to report per-variant purity and ploidy. However, if you target a certain stage specifically and exclude others, umccrise will use the target stage without the dependency. E.g. `germline` will not contain germline leakage from somatic variants if `somatic` stage is not included, `coverage` report will not pull `purple` and will use reasonable purity threshold defaults instead, `multiqc` will report on available QC data only, etc.

- Support reference data location (or prefix) in form of an `s3://` or `gds://` URL. E.g.

```
umccrise /input --genomes s3://umccr-refdata-dev/genomes
```

Would check `s3://umccr-refdata-dev/genomes_102`, `s3://umccr-refdata-dev/genomes_10`, and `s3://umccr-refdata-dev/genomes`, assuming that the [reference_data](https://github.com/umccr/reference_data) package version is `1.0.2`.

Downloads the reference data locally into a `~/umccrise_genomes`.

--------------------

0.17.12, 0.17.11 (4 May 2020)

- Dropping sv-prioritize filters from bcbio, rescuing SVs wrongly marked as "intergenic"
- CPSR: hiding the genes pane which causes random bugs

0.17.10, 0.17.9 (26 Mar 2020)

- Fix SV processing for FFPE

0.17.8 (23 Mar 2020)

- Fix the core usage inside of a Docker image, so COBALT and AMBER will benefit from multiple cores and AWS runs will get much faster now (thanks Peter for spotting the issue)

0.17.7 (14 Mar 2020)

- CPSR: when prioritizing a mutation effects, always pick the transcripts where the mutations affects the protein.

0.17.6 (12 Mar 2020)

- To make sure CodeBuild builds code from branches, modify install.sh so it
doesn't clone the repo, and move Dockerfile to the root to enable copying
the source.

0.17.5 (12 Mar 2020)

- Fix germline variants count in MultiQC

0.17.4 (11 Mar 2020)

- Pin PCGR and CPSR and show their correct versions.
- Pin MultiQC to vladsaveliev/MultiQC rather than umccr-illumina/MultiQC

0.17.3 (9 Mar 2020)

- Fix https://github.com/umccr/umccrise/issues/41

0.17.2 (6 Mar 2020)

Fixes
- Fix PURPLE memory requirements for single-core runs
- Remove entrypoint from the Dockerfile

0.17.1 (27 Feb 2020)

- Fix PCGR and CPSR transcript prioritization issues https://github.com/sigven/cpsr/issues/26
- Output SAGE results into {batch}/small_variants/sage/{batch}-sage.vcf.gz

0.17.0 (1 Feb 202)

- Massive cancer report updates by Peter Diakumis:
    - Added new somatic signatures
    - SV tables restrucurised and got much more readable
- Support DRAGEN project on input instead of bcbio

--------------------

0.16.0 (31 Dec 2019)

- Add CodeBuild CI support (thanks @brainstorm for the great work!)
- Add hg38 MultiQC background samples
- Update PURPLE, add new plots/stats into the cancer report
- FFPE pipeline: remove purple-inferred SVs
- SV prioritization: move filtering BNDs with SR > PR from Rmd into the prioritization script (so this filter reflects in VCF)
- Added GRIDSS/PURPLE/LINX wrapper (`gpl`) and env (`umccrise_hmf`)
- hpc_utils: refactoring of genomes_dir detection functions. Add set_genomes_dir
- Fix multithreading when running in a Docker
- Support Gadi
- Update PCGR to 0.8.4 and CPSR to 0.5.2:
    - VEP updated to 98.3
    - Data bundle updates (CIViC, ClinVar, CancerMine, UniProt)
    - Fixed multiple bugs in hg38 support for CPSR:
        - Fixed transcript priority issue (https://github.com/sigven/pcgr/issues/79)
        - Fixed Germline biomarker parsing bugs (https://github.com/sigven/pcgr/pull/85, https://github.com/sigven/pcgr/commit/e4e428eafa8f2a93fa7053908a447ef064d5a9c9)
        - Fixed https://github.com/sigven/cpsr/issues/15 with VEP update
        - LofTee plugin fixed
    - Conda build optimizations: added LofTee and perl-bio-big-file (https://github.com/sigven/pcgr/pull/98)
    - Fixed crosstalk and BSgenome issues https://github.com/sigven/cpsr/issues/19
- Support unprioritized SV VCFs naming (when sv-prioritize is disabled in bcbio)

--------------------

0.15.6 (28 Aug 2019)

- SV prioritization: 
    - prioritize CN=0 events in tumor suppressors genes
    - fix transcript order in PCGR (use Rank first)
    - report all (principal+alt) transcripts comma-separated
    - skip interaction events
- Cancer report SV table: 
    - Add SV event number
    - Fix ref/alt read support labels
- Cancer report: reorganize blocks, increase report width
- SV VCF: split samples, keep tumor (for Pierian)
- Add MultiQC Docker container for bcl2fastq and interop

0.15 (21 May 2019)

- Update PCGR to 0.8.1 and CPSR to 0.4.1 (with a new scheme for germline filtering)
- Pass tumor purity and ploidy into PCGR
- VEP annotation pick order: pick APPRIS first instead of canonical (to match with SV)
- Add the detailed [workflow description](https://github.com/umccr/umccrise/blob/master/workflow.md)

--------------------

0.14 (1 Apr 2019)

- Update SV prioritization method
- New SV table in the rmd report
- New PURPLE tools, PURPLE annotates and rescues SVs
- Panel of normals rebuilt with lowfreq tumor-only Mutect2, threshold is raised to 5
- Coverage report threshold: adjust to purity (CACAO and mosdepth)
- Support of Cromwell MultiQC
- Report somatic germine leakage (gnomAD variants) in CPSR

--------------------

0.13 (15 Feb 2019)

- Re-enable [CACAO](https://github.com/vladsaveliev/cacao)
- Support Cromwell runs
- MultiQC:
  - add somatic and germline filtering stats
  - run bcftools-stats for filtered calls
  - sort samples in tumor-normal order
  - add columns: secondary alignments, MQ=0 reads, viral
  - clean up
- Integrate [SAGE](https://github.com/hartwigmedical/hmftools/tree/master/sage): to refine variant calls in hotspot locations
- SVs:
  - Prioritize with [simple_sv_annotation](https://github.com/AstraZeneca-NGS/simple_sv_annotation) 
    using [HMF fusions](https://nc.hartwigmedicalfoundation.nl/index.php/s/a8lgLsUrZI5gndd?path=%2FHMF-Pipeline-Resources) 
    and UMCCR cancer genes
- FFPE:
  - Hard subset small variants to non-gnomad if >500k, subset to cancer genes if doesn't help
  - Hard subset SVs to non-intergenic if >1000, subset to prioritized if doesn't help
  - Rerun BPI on priority SVs when disabled in bcbio
- Cancer report significantly updated (thanks to Peter's work!):
  - Add mutational signatures descriptions
  - Discontinue bookdown for a speed boost
  - Restructure tables and plots
  - Add explanations how files are generated in the pipeline

--------------------

0.12 (23 Jan 2019)

- Somatic annotation/filtering:
    - Filter VarDict strand bias and HP indels
    - Abstract somatic annotation into somatic_vcf_annotate script
    - Add more sources: HMF mappability, ICGC, COSMIC, TCGA hotspots, ENCODE, LCR
- Populate VD and AD in final VCFs (update to pcgr_prep)
- Support Strelka2-only runs, populate Strelka2 AD for Purple
- Correct threading for cluster runs
- Key genes:
    - Update umccr cancer genes to 2019 (remove TSG, add cancermine and NGC)
    - Use NGC+cancermine oncogene/tsg defenitions
    - Add predispose BED files for germline runs
- Update Purple to sambamba-free Amber
- Integrate [Conpair](https://github.com/vladsaveliev/Conpair) (replaces VerifyBAMID and Peddy)

--------------------

0.11 (6 Dec 2018)

- Add Per-sample MultiQC reports with reference samples
- Add Sig.rmd Purple table
- Add summary Purple-vs-bcbio comparison

--------------------

0.10 (30 Nov 2018)

- SNV filtering improvements: 
    - Assuming calling with 1%
    - First round of PoN filter: remove with >=3 hits
    - Annotate with PCGR and keep all tier 1-3
    - Tier 4 and noncoding - filter with: 
        - gnomad_AF <=0.02
        - PoN hits >=1
        - indels in bad promoters
        - DP<25 & AF<5% in tricky regions or GiaB LCR
- Integrate CACAO (temporary disabled)

--------------------

0.9 (15 Nov 2018)

- Integrate Purple
- Integrate CPSR (PCGR reports for germline variants)
- Update PCGR
- Use Arthur's [new list](https://trello.com/c/7j3KFMiL/184-umccr-cancer-genes) of cancer key genes for minibams
    and coverage
- Use Sigve's+Sean's predispodition gene lists for germline reporting

--------------------

0.8 (7 Aug 2018)

- SV filtering: require all of FILTER fields to be PASS or rejected by sv-prioritizer
- SV filtering: always keep SR or PR above 10x regardless of BPI_AF (but above 5x only if BPI_AF>10%)
- SV filtering: check BPI_AF in both breakends
- Add option --no-s3 to skip IGV upload; use it in tests
- Restart failed jobs tree times (helps with PCGR random issues)
- Exlcuding/including samples: perform before check for correct YAML
- Exlcuding/including samples: support {batch}__{sample} identifiers
- Sopport the new and old germline VCF bcbio naming (<normal>-<caller>.vcf.gz and <normal>-germline-<caller>.vcf.gz)

--------------------

0.7 (26 Jun 2018)

- Dockerized version with --docker
- Custom reference data with --pon, --ref-fasta, --truth-regions, or --bcbio-genomes

--------------------

0.6 (15 May 2018)

- Use snakemake groups to optimize execution on cluster
- Use submission wrapper to automate cluster resources
- Filter structural variants: BPI_AF>=10% and read support>=5
- Keep sv-prioritize hard-filtered variants
- 0.6.1
  - Add circos plots into the Rmd report. Generated by https://github.com/umccr/rock

--------------------

0.5 (7 May 2018)

- PCGR is deployed directly on Spartan, so no AWS dependency.
- Add PCGR wrapper: `pcgr variants.vcf.gz cnv.tsv -o results [-g hg38]`
- Correctly providing memory resources on HPC to avoid oom-kill
- On Spartan, support `--cluster-auto` to automatically substitute proper cluster parameters

--------------------

0.4 (17 Apr 2018)

- Propagate snakemake's cluster options to the wrapper
- Propagate snakemake's `--rerun-incomplete` and to the wrapper
- PCGR: stop using the UID, skip uploading if already uploaded
- Fix `HTTP_PROXY` setting for running on Spartan worker nodes
- New panel of normals: add the large A5 cohort, normalize normals, check only POS for indels
- Correctly set threads to mibibams samtools runs

--------------------

0.3 (17 Jan 2018)

- Refactor output folder structure (see docs)
- `pcgr_download` taget to automatically pull results
- Automatically attempt downloading PCGR results after minibams
- Minibams always run in the ned
- Use `--unlock` to automate the snakemake unlock pattern to continue unterrupted runs
- Move `vcf_stuff` to a separate repo
- Remove GL* chromosomes from the cnvkit diagram
- PCGR tomls: changed mutsignatures_normalization to genome
- Rmd: add SV-prioritize table
