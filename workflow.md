# Umccrise workflow


## SNPs and small indels

### Somatic

Umccrise post-processes somatic variants in order to remove most of the artefacts and germline leakage, but at the same time be conservative enough in known hotspot loci to make sure to preserve actionable variants. The following post-processing scheme is designed for that purpose:

1. Take passing "ensemble" somatic VCF from [bcbio](https://github.com/umccr/workflows/tree/master/bcbio). "Ensemble" has variants supported by at least 2 of 3 callers (we use [strelka2](https://github.com/Illumina/strelka), [vardict](https://github.com/AstraZeneca-NGS/VarDict), and [mutect2](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.4/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php)). Unless it's an FFPE sample, in which case usually only strelka2 is used. Variants are called at AF>1%.
2. Sort VCF by coordinate, extract `PASS` calls (which applies only to strelka2, as ensemble is already `PASS`-only).
3. Run [SAGE](https://github.com/hartwigmedical/hmftools/tree/master/sage) and add the result to the VCF. SAGE is a low-frequency variant caller with a high precision, created by Hartwig Medical Foundation. Instead of the whole genome, it targets only coding regions for inframe indels, and [known hotspot sites](https://raw.githubusercontent.com/umccr/workflows/master/genes/hotspots/hotspots.tsv) from the following list:
	*  [Cancer Genome Interpreter](https://www.cancergenomeinterpreter.org/home) 
	*  [CIViC](http://civic.genome.wustl.edu/) - Clinical interpretations of variants in cancer
	*  [OncoKB](https://oncokb.org/) - Precision Oncology Knowledge Base
4. Annotate the VCF against the reference sources:
	*  SAGE hotspots (CGI, CiViC, OncoKB), 
	*  GiaB confidence regions,
	*  [GnomAD](https://gnomad.broadinstitute.org) whole genome common variants (max population frequency > 1%),
	*  Low complexity regions (Heng Li's)
	*  LCR, low and high GC regions, self-chain and bad promoter regions (GA4GH),
	*  [ENCODE blacklist](https://github.com/Boyle-Lab/Blacklist),
	*  Segmental duplication regions (UCSC),
	*  UMCCR [panel of normals](https://github.com/umccr/vcf_stuff/blob/master/vcf_stuff/panel_of_normals/story/panel_of_normals.md), build from tumor-only mutect2 calls from ~200 normal samples.

[ ] _TODO: link to PanelApp or original data where possible._

5. If after removing non-hotspot GnomAD variants there are still > 500k somatic variants left flag the sample as highly mutated (or FFPE) and limit all calls to to cancer genes only (to avoid downstream PCGR and UMCCR cancer reports choking).

[ ] _TODO: referfence "cancer genes" with a link to the source or PanelApp._

6. Standardize the VCF fields: add new `INFO` fields `TUMOR_AF`, `NORMAL_AF`, `TUMOR_DP`, `NORMAL_DP`, `TUMOR_VD`, `NORMAL_VD` (for use with PCGR), and `AD FORMAT` field (for use with PURPLE).
7. Run [PCGR](https://github.com/sigven/pcgr) to annotate VCF with more external sources:
	* [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html)  to infer the functional impact
	* TCGA and ICGC recurrent variants
	* [Open Targets Platform](https://targetvalidation.org/)
    * [ClinVar](http://www.ncbi.nlm.nih.gov/clinvar/) - Database of variants with clinical significance
    * [CancerMine](https://zenodo.org/record/2662509#.XM0xMdMzaL4) - Literature-derived database of tumor suppressor genes/proto-oncogenes
    * [DoCM](http://docm.genome.wustl.edu) - Database of curated mutations
    * [CBMDB](http://www.cancergenomeinterpreter.org/biomarkers) - Cancer Biomarkers database 
    * [DisGeNET](http://www.disgenet.org) - Database of gene-tumor type associations
    * [Cancer Hotspots](http://cancerhotspots.org) - Resource for statistically significant mutations in cancer
    * [dBNSFP](https://sites.google.com/site/jpopgen/dbNSFP) - Database of non-synonymous functional predictions
    * [UniProt/SwissProt KnowledgeBase](http://www.uniprot.org) - Resource on protein sequence and functional information
    * [Pfam](http://pfam.xfam.org) - Database of protein families and domains
    * [DGIdb](http://dgidb.genome.wustl.edu) - Database of targeted cancer drugs
    * [ChEMBL](https://www.ebi.ac.uk/chembl/) - Manually curated database of bioactive molecules
   
    PCGR also classifies by tiers based on annotations and functional impact.  At the end, this step adds `INFO` fields into the VCF: `TIER`, `SYMBOL`, `CONSEQUENCE`, `MUTATION_HOTSPOT`, `INTOGEN_DRIVER_MUT`, `TCGA_PANCANCER_COUNT`, `CLINVAR_CLNSIG`, `ICGC_PCAWG_HITS`, `COSMIC_CNT`.

[ ] _TODO: Harmonize gene lists with overall workflow description [gene lists](https://github.com/umccr/workflows#gene-lists)._

8. [Filter variants](https://github.com/umccr/vcf_stuff/blob/master/scripts/filter_somatic_vcf) to remove germline variants and artefacts, but rescue known hotspots/actionable variants:
	*  Rescue SAGE hotspots (CGI, CiViC, OncoKB),
	*  Rescue PCGR TIER 1 and 2 (strong and potential clinical significance, according to [ACMG](https://www.ncbi.nlm.nih.gov/pubmed/27993330) standard guidelines,
	*  Additionally rescue all driver mutations ([Intogen](https://www.intogen.org/)); [mutation hotspots](http://cancerhotspots.org/]); [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) likely pathogenic or uncertain; COSMIC hits >=10; TCGA pancancer hits>=5; ICGC PCAWG hits >= 3 (all annotated by PCGR).
	*  For all non-rescued variants, apply LCR, PoN, depth and AF filters. Remove variants for which one or more of the following conditions apply:
		* `AF`<10%,
		* Common variant in GnomAD (max `population AF`>=1%), add into the germline set (see below),
		* Present in >=5 samples of the Panel of Normal set,,
		* InDel in a "bad promoter" regions,
		* Remove variants overlapping the ENCODE blacklist,
		* `VD`<4 (or <6 in tricky regions: GC<15, GC>70, LCR (GA4GA and Heng's Li), GiaB non-confident, segmental duplications), 
		* VarDict strand biased variants (single strand support for ALT, while REF has both; or REF and ALT have opposite supporting strands), unless supported by all other callers,
		* Dynamic AF threshold for indels in MSI for VarDict (unless supported by all other callers).

[ ] _TODO: Define what "rescue" means in this context (set to non-filtered?). Should the order be different? First filter, then rescue anything that was filtered but is present in a hotspot or is PCGR Tier1/2? Can we define what a 'hit' in COSMIC, TCGA, etc., means? Define "bad promoter" (again, if needed). The "dynamic AF" threshold might need additional information. The "tricky" region file needs to be linked (BED), it differs slightly from the list defined above which does not include GC extremes._

9. Report passing variants using [PCGR](https://github.com/sigven/pcgr), classified by the ACMG tier system.


### Germline

The idea is to simply bring germline variants in cancer predisposition genes:

1. Take passing "ensemble" germline VCF from [bcbio](https://github.com/umccr/workflows/tree/master/bcbio). "Ensemble" has variants supported by at least 2 of 3 callers (we use strelka2, vardict, and GATK Haplotype Caller).
2. Sort VCF by coordinate, extract `PASS` calls.
3. Add back variants from somatic calling filtered as common GnomAD.
4. Subset variants to a list of ~200 cancer predisposition genes, which is build by [CPSR](https://github.com/sigven/cpsr) from 3 curated sources: [TCGA](https://www.ncbi.nlm.nih.gov/pubmed/29625052) pan-cancer study, [COSMIC CGC](https://cancer.sanger.ac.uk/census), and [Norwegian Cancer Genomics Consortium](http://cancergenomics.no/).
5. Report variants using CPSR, which classifies variants in the context of cancer predisposition by overlapping with [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) pathogenic and VUS variants and GnomAD rare variants. It also ranks variants according pathogenicity score by ACMG and cancer-specific criteria.

[ ] _TODO: Confirm the changes in here, please. Does this still match the updates to preferred transcript choice? Do we want to mention the logic within CPSR on how it picks the right transcript? Again, link to workflow document for gene lists, or harmonize._

[ ] _TODO: Old text from Google Docs below. Is this still valid? If so, merge in:_

> * ClinVar variants: pre-classified variants according to a five-level tier scheme (Pathogenic to Benign):
> * Non-ClinVar variants: classified by CPSR according to a five-level tier scheme (Pathogenic to Benign)
> * Secondary findings (optional) - pathogenic ClinVar variants in the ACMG recommended list for reporting of incidental findings,
> * GWAS hits (optional) - variants overlapping with previously identified hits in genome-wide association studies (GWAS) of cancer phenotypes (i.e. low to moderate risk conferring alleles), using NHGRI-EBI Catalog of published genome-wide association studies as the underlying source.
> 
> The unclassified non-ClinVar variants are assigned a pathogenicity level based on the aggregation of scores according to previously established ACMG criteria. The ACMG criteria includes cancer-specific criteria, as outlined and specified in several previous studies.

## Structural variants

The idea is to report gene fusions, exon deletions, high impact and LoF events in tumor suppressors, and prioritize events in cancer genes.

1. Start with the somatic SV VCF from [bcbio](https://github.com/umccr/workflows/tree/master/bcbio) called by [Manta](https://github.com/illumina/manta) SV caller. 
2. Refine SVs using Hartwig break-point-inspector, which locally re-assembles SV loci to get more accurate breakpoint positions and AF estimates.
3. Filter low-quality calls:
   * require split or paired reads support at least 5x,
   * for low frequency variants (<10% at both breakpoints), require read support 10x,
   * require paired reads support to be higher than split read support for BND events
4. Annotate variants impact using [SnpEff](http://snpeff.sourceforge.net/SnpEff_manual.html) according to the Ensembl gene model and [Sequence ontology](http://www.sequenceontology.org) terminology.
5. Subset annotations to [APPRIS principal transcripts](http://appris.bioinfo.cnio.es/#/), keeping one main isoform per gene.
6. Use variants as a guidance for PURPLE CNV calling (see below). PURPLE will adjust and recover breakpoints at copy number transitions, and adjust AF based on copy number, purity and ploidy.
7. Prioritize variants with [simple_sv_annotation](https://github.com/vladsaveliev/simple_sv_annotation) on a 4 tier system - 1 (high) - 2 (moderate) - 3 (low) - 4 (no interest):
    * exon loss
       * on cancer gene list (1)
       * other (2)
    * gene fusion
       * paired (hits two genes)
          * on list of known pairs (1) (curated by [HMF](https://resources.hartwigmedicalfoundation.nl))
          * one gene is a known promiscuous fusion gene (1) (curated by [HMF](https://resources.hartwigmedicalfoundation.nl))
          * on list of [FusionCatcher](https://github.com/ndaniel/fusioncatcher/blob/master/bin/generate_known.py) known pairs (2)
          * other:
             * one or two genes on cancer gene list (2)
             * neither gene on cancer gene list (3)
       * unpaired (hits one gene)
           * on cancer gene list (2)
           * others (3)
    * upstream or downstream (a specific type of fusion, e.g. one gene is got into control of another gene's promoter and get over-expressed (oncogene) or underexpressed (tsgene))
       * on cancer gene list genes (2)
    * LoF or HIGH impact in a tumor suppressor
       * on cancer gene list (2)
       * other TS gene (3)
    * other (4)    

[ ] _TODO: Compare against this version: which one is correct / clearer?_

> a. If the feature type is "interaction" (protein-protein interactions), or a generic "sequence_feature", discard the event
> b. If there is an associated transcript ID and it's not in list of principal or alternative APPRIS transcripts, discard the event
> c. exon loss:
>   i. of a umccr gene (tier 2)
>   ii. other (3)
> d. gene fusion:
>   i. paired (hits exactly two genes)
>      * on list of known pairs (tier 1) (curated by [HMF](https://resources.hartwigmedicalfoundation.nl))
>      * one gene is a known promiscuous fusion gene (1) (curated by [HMF](https://resources.hartwigmedicalfoundation.nl))
>      * on list of FusionCatcher known pairs (2)
>      * one or two genes are on umccr gene list (2)
>   ii. unpaired (hits one gene)
>      * on umccr gene list (2)
> e. upstream or downstream: same logic as with fusions, but 1 tier down (i.e. tier 2->3, 1->2)
> f. HIGH impact, and hits a tumor suppressor or an entire chromosome
>   i. chromosome (2)
>   ii. key tumor suppressor (2)
>   iii. other tumor suppressors (3)
>   iv. if the copy number is zero, upgrade 1 tier up (i.e. 2->1)
> g. any other event on umccr gene list (3)
> h. all other events (4)

8. If the number of events is over 10k (e.g. FFPE), keep only tiers 1-2-3
9. For tiers 3-4:
  a. require split or paired reads support of at least 5x,
  b. for low frequency variants (<10% at both breakpoints), require read support of at least 10x,
  c. require paired reads support to be higher than split read support for BND events
10. Unless FFPE, feed the variants into PURPLE to recover SVs from copy number transitions
11.  Report tiered variants in the UMCCR cancer report.

[ ] _TODO: BPI is the problematic bit here. Have to rely on us shipping it / can't link to a source. Maybe just cite the Hartwig paper. Peter, can you go over this section and compare against your detailed notes? Same as with other sections: harmonize gene lists and link to actual data (e.g., APPRIS transcripts)_


## Somatic CNV

The idea is to report significant CNV changes in key cancer genes and disruptions in tumor suppressors. And also calculate sample purity and ploidy profile.

We almost entirely rely on Hartwig's [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) workflow in this step. The PURPLE pipeline outlines as follows:

* Calculate B-allele frequencies (BAF) using AMBER subworkflow,
* Calculate read depth ratios using COBALT subworkflow,
* Perform segmentation (uses structural variant breakpoints for better guidance),
* Estimate the purity and copy number profile (uses somatic variants for better fitting), 
* Plot a circos plot that visualizes the CN/ploidy profile, as well as somatic variants, SVs, and BAFs,
* Rescue structural variants in copy number transitions and filter single breakends,
* Estimate overall tumor samples purity range,
* Determine gender, 
* Report QC status of the sample, that will fail if the structural variants do not correspond to CN transitions, and gender is inconsistently called from BAFs and from the coverage.

From the PURPLE output, we report in the cancer report:

* Circos plot
* Minimal and maximal copy numbers in key cancer genes, that indicate amplifications/deletions as well as CN transitions that should match SVs
* QC status
* We also use Purity to adjust coverage reporting thresholds.
* Genome-wide CNV segments with breakpoint information. Includes segment CN, minor/major allele ploidy, type of SV support at start/end of segment, CN determination method, BAF/BAF count, GC%, Cobalt window coverage


## MultiQC

MultiQC aggregates QC from different tools. We report the following:
	
* Sample contamination level (for both tumor and normal) and tumor/normal concordance (by [Conpair](https://github.com/nygenome/Conpair)),
* Ancestry and sex (by [Peddy](https://github.com/brentp/peddy)),
* Mapping QC: the number of mapped reads, paired reads, secondary or duplicated alignments, average coverage (using samtools stats and mosdepth in bcbio),
* Viral DNA content (in bcbio),
* Number of pre- and post-filtered SNPs and indels (by Umccrise) which indicates germline leakage,
* Coverage profile by goleft,
* Variants QC for filtered germline and somatic variants (by bcftools),
* Reads QC (by FastQC).

We also include reference "good" samples in the background for comparison.


## Coverage

The idea is to see if we can miss variants due to abnormal coverage (e.g. because of copy numbers or abnormal GC).

* Run [`goleft indexcov`](https://github.com/brentp/goleft/tree/master/indexcov) to plot fast global coverage overview from a BAM or CRAM index.
* Run [`mosdepth`](https://github.com/brentp/mosdepth) to calculate quantized coverage in exons of cancer genes if interest, using 4 groups for quantization: NO_COVERAGE, LOW_COVERAGE, CALLABLE, HIGH_COVERAGE. For tumor, the thresholds are adjusted by average purity (as reported by PURPLE). The low coverage threshold is 12x divided by purity (the minimal coverage to call a pure heterozugous variant), the high coverage threshold is 100 divided by purity (suspiciously high coverage, might indicate an LCR, a repeat, or a CN).
* Run [CACAO](https://github.com/sigven/cacao) using the same thresholds to calculate coverage in loci of interest. For germline variants, use [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar) predisposition variants. For somatic variants, use [CiViC](https://civicdb.org/) actionable variants and [cancerhotspots.org](https://www.cancerhotspots.org/) somatic hotspots. This step generates two HTML reports (one for somatic, one for germline variants).


## Reports

Umccrise produces 6 reports:

* [PCGR](https://github.com/sigven/pcgr) containing small somatic variants (SNPs and indels) classified according to [ACMG](https://www.ncbi.nlm.nih.gov/pubmed/27993330) guidelines, and MSI status of the sample.
* [CPSR](https://github.com/sigven/cpsr) containing small germline variants (SNPs and indels) in cancer predisposition genes, ranked by ACMG guidelines and cancer-specific criteria.
* [CACAO](https://github.com/sigven/cacao) for tumor sample, reporting coverage for clinically actionable and pathogenic loci in cancer
* CACAO for normal sample, reporting coverage in likely pathogenic variants cancer predisposition protein-coding genes
* [MultiQC](https://multiqc.info) report with QC stats and plots
* UMCCR cancer report containing:
	* Somatic mutation profile (global and in cancer genes),
	* Mutational signatures (by the MutationalPatterns R package),
	* Structural variants,
	* Copy number variants,
	* PURPLE QC status,
	* Circos plot.


## Minibam

We generate a portable mini-version of the BAM file by subsetting it to key cancer genes, somatic variants, and SVs breakpoints.


## Key cancer genes

For reporting and variant prioritization, we prepared a [UMCCR cancer key genes set](https://github.com/vladsaveliev/NGS_Utils/blob/master/ngs_utils/reference_data/key_genes/make_umccr_cancer_genes.Rmd). It's build of off several sources:

* Cancermine with at least 2 publication with at least 3 citations,
* NCG known cancer genes,
* Tier 1 COSMIC Cancer Gene Census (CGC),
* CACAO hotspot genes (curated from ClinVar, CiViC, cancerhotspots),
* At least 2 matches in the following 5 sources and 8 clinical panels:
   * Cancer predisposition genes (CPSR list),
   * COSMIC Cancer Gene Census (tier 2),
   * AZ300, 
   * Familial Cancer, 
   * OncoKB annotated,
   * MSKC-IMPACT, 
   * MSKC-Heme, 
   * PMCC-CCP, 
   * Illumina-TS500, 
   * TEMPUS, 
   * Foundation One, 
   * Foundation Heme, 
   * Vogelstein.

The result is a list of 1248 genes.

[ ] _TODO: Harmonize gene list._
