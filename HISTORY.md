0.16 (TBR)

- Add CodeBuild CI support (thanks @brainstorm for the great work!)
- Add hg38 MultiQC background samples
- Update PURPLE, add new plots/stats into the cancer report
- FFPE pipeline: remove purple-inferred SVs
- SV prioritization: move filtering BNDs with SR > PR from Rmd into the prioritization script (so this filter reflects in VCF)
- Added GRIDSS/PURPLE/LINX wrapper (`gpl`) and env (`umccrise_hmf`)
- hpc_utils: refactoring of genomes_dir detection functions. Add set_genomes_dir

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

0.14 (1 Apr 2019)

- Update SV prioritization method
- New SV table in the rmd report
- New PURPLE tools, PURPLE annotates and rescues SVs
- Panel of normals rebuilt with lowfreq tumor-only Mutect2, threshold is raised to 5
- Coverage report threshold: adjust to purity (CACAO and mosdepth)
- Support of Cromwell MultiQC
- Report somatic germine leakage (gnomAD variants) in CPSR

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

0.11 (6 Dec 2018)

- Add Per-sample MultiQC reports with reference samples
- Add Sig.rmd Purple table
- Add summary Purple-vs-bcbio comparison

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

0.9 (15 Nov 2018)

- Integrate Purple
- Integrate CPSR (PCGR reports for germline variants)
- Update PCGR
- Use Arthur's [new list](https://trello.com/c/7j3KFMiL/184-umccr-cancer-genes) of cancer key genes for minibams
    and coverage
- Use Sigve's+Sean's predispodition gene lists for germline reporting

0.8 (7 Aug 2018)

- SV filtering: require all of FILTER fields to be PASS or rejected by sv-prioritizer
- SV filtering: always keep SR or PR above 10x regardless of BPI_AF (but above 5x only if BPI_AF>10%)
- SV filtering: check BPI_AF in both breakends
- Add option --no-s3 to skip IGV upload; use it in tests
- Restart failed jobs tree times (helps with PCGR random issues)
- Exlcuding/including samples: perform before check for correct YAML
- Exlcuding/including samples: support {batch}__{sample} identifiers
- Sopport the new and old germline VCF bcbio naming (<normal>-<caller>.vcf.gz and <normal>-germline-<caller>.vcf.gz)

0.7 (26 Jun 2018)

- Dockerized version with --docker
- Custom reference data with --pon, --ref-fasta, --truth-regions, or --bcbio-genomes

0.6 (15 May 2018)

- Use snakemake groups to optimize execution on cluster
- Use submission wrapper to automate cluster resources
- Filter structural variants: BPI_AF>=10% and read support>=5
- Keep sv-prioritize hard-filtered variants
- 0.6.1
  - Add circos plots into the Rmd report. Generated by https://github.com/umccr/rock

0.5 (7 May 2018)

- PCGR is deployed directly on Spartan, so no AWS dependency.
- Add PCGR wrapper: `pcgr variants.vcf.gz cnv.tsv -o results [-g hg38]`
- Correctly providing memory resources on HPC to avoid oom-kill
- On Spartan, support `--cluster-auto` to automatically substitute proper cluster parameters

0.4 (17 Apr 2018)

- Propagate snakemake's cluster options to the wrapper
- Propagate snakemake's `--rerun-incomplete` and to the wrapper
- PCGR: stop using the UID, skip uploading if already uploaded
- Fix `HTTP_PROXY` setting for running on Spartan worker nodes
- New panel of normals: add the large A5 cohort, normalize normals, check only POS for indels
- Correctly set threads to mibibams samtools runs

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
