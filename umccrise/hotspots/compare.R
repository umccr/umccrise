########### Compare hotspots
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(purrr)

(hart_hs = read_tsv('/Users/vsaveliev/bio/genomes/GRCh37/hotspots/KnownHotspots.tsv', 
                   col_names = c('CHROM', 'POS', 'REF', 'ALT'),
                   col_types = 'cicc') %>% 
  mutate(source = 'hart'))

(pcgr_hs = read_tsv('/Users/vsaveliev/bio/genomes/pcgr/data/grch37/cancer_hotspots/cancer_hotspots.tsv',
                    col_types = 'ciccc') %>% 
    mutate(source = 'pcgr') %>% 
    separate(INFO, c("change", "transcript", "cancertype"), sep=";") %>% 
    mutate(
      change = str_replace(change, "MUTATION_HOTSPOT=", ""),
      transcript = str_replace(transcript, "MUTATION_HOTSPOT_TRANSCRIPT=", ""),  # TODO: filter with canonical transcripts
      cancertype = str_replace(cancertype, "MUTATION_HOTSPOT_CANCERTYPE=", "")
    ) %>% 
    select(-transcript)
)


################ Comparing HMF and PCGR ################

# A lot of repeated SNPs
# full_join(hart_hs, pcgr_hs, by=c("CHROM", "POS", "REF", "source")) %>% 
#   group_by(CHROM, POS, REF, source) %>% 
#   count(sort = T)

# Merging finally
full_join(
  hart_hs %>% mutate(INFO = 'HMF') %>% select(-source), 
  pcgr_hs %>% mutate(INFO = str_c(change, cancertype, sep = "|")) %>% select(-change, -cancertype, -source), 
  by=c("CHROM", "POS", "REF", "ALT", "INFO")
) %>% 
  group_by(CHROM, POS, REF, ALT) %>% 
  summarise(INFO = str_c(INFO, collapse = ";")) %>% 
  print(n = 100)

library(nVennR)
plotVenn(list(hart = pcgr_hs %>% select(POS), pcgr = pcgr_hs %>% select(POS)))


(overlap <- intersect(hart_hs %>% select(CHROM, POS, REF, ALT), 
                      pcgr_hs %>% select(CHROM, POS, REF, ALT)))

pcgr_hs %>% count()
hart_hs %>% count()
(pcgr_locs = pcgr_hs %>% distinct(CHROM, POS))
(hart_locs = hart_hs %>% distinct(CHROM, POS))

pcgr_locs %>% count()
hart_locs %>% count()
intersect(hart_locs, pcgr_locs) %>% count()
setdiff(hmf_locs, pcgr_locs) %>% count()
setdiff(pcgr_locs, hmf_locs) %>% count()



################ Adding CACAO ################

cacao_dir = "/Users/vsaveliev/git/cacao/data/"
sources = c(
  "civic" = "cacao.civic.grch37.tsv",
  "clinvar_path" = "cacao.clinvar_path.grch37.tsv",
  "hotspot" = "cacao.hotspot.grch37.tsv"
)
(cacao_hs_comb = sources %>%
    map2(names(sources), (function(fp, name) { 
      fp %>% 
        str_c(cacao_dir, .) %>%
        read_tsv(col_types = cols(chrom = "c")) %>% 
        select(chrom, start, end) %>%
        mutate(source = name)
    })) %>% 
    bind_rows()
)

(cacao_hs_all = read_tsv(str_c(cacao_dir, "cacao.grch37.bed"), col_types = "ciic", col_names = c("chrom", "start", "end", "ann")) %>% 
    mutate(
      source = str_split(ann, ":") %>% map_chr(~ .[[1]]),
      source = factor(source, levels = c("civic", "clinvar_path", "hotspot"))
    )
)

cacao_hs_comb %>% count(source)
# civic         1095
# clinvar_path 16801
# hotspot       1024

cacao_hs_all %>% count(source)
# civic         1092
# clinvar_path 16801
# hotspot       1024
# NA               3

# Sources mostly match, but 3 variants (really, 2) in the all file missing a source, and it looks like it's civic:
cacao_hs_all %>% filter(is.na(source))
# 5      1295222  1295233 TERT:upstream_gene_variant_C228T        NA    
# 5      1295222  1295233 TERT:upstream_gene_variant_C228T        NA    
# 9     21975011 21975022 CDKN2A:upstream_gene_variant_CA12951936 NA    

# Checking:
cacao_hs_comb %>% inner_join(cacao_hs_all %>% filter(is.na(source)), by = c("start"))

# Yep. So just fixing the BED file.
# /Users/vsaveliev/git/cacao/data$ sort -k1,1V -k2,2n cacao.grch37.bed | uniq -c > cacao.grch37.bed.unique                                                                                                                                                                                                                                                                                                                                                                                                                                   *[master]
# /Users/vsaveliev/git/cacao/data$ mv cacao.grch37.bed.unique cacao.grch37.bed
# Same for hg19 and hg38


################ Comparing HMF, PCGR and CACAO ################

# We want to feed SAGE with the same hotspots as PCGR and CACAO. But we also want to report coverage of HMF hotspots and prioritize them in PCGR.
# For that, we need to merge all the lists.
# The problem is that they have different representations. PCGR and HMF use POS-REF-ALT, but CACAO uses START-END.
# Ask Sigve?



# Now comparing coordinates with HMF and PCGR hotspots
full_join(
  cacao_hs_all %>% rename(CHROM = chrom, POS = start, END = end, INFO = source) %>% select(-ann),
  hart_hs %>% mutate(INFO = 'HMF') %>% select(-source), 
  pcgr_hs %>% mutate(INFO = str_c(change, cancertype, sep = "|")) %>% select(-change, -cancertype, -source), 
  by=c("CHROM", "POS")
) %>% 
  group_by(CHROM, POS) %>% 
  # summarise(INFO = str_c(INFO, collapse = ";")) %>% 
  summarize(n = n()) %>% 
  arrange(desc(n)) %>% 
  print(n = 100)

library(nVennR)
plotVenn(list(hart = pcgr_hs %>% select(POS), pcgr = pcgr_hs %>% select(POS)))


(overlap <- intersect(hart_hs %>% select(CHROM, POS, REF, ALT), 
                      pcgr_hs %>% select(CHROM, POS, REF, ALT)))

pcgr_hs %>% count()
hart_hs %>% count()
(pcgr_locs = pcgr_hs %>% distinct(CHROM, POS))
(hart_locs = hart_hs %>% distinct(CHROM, POS))

pcgr_locs %>% count()
hart_locs %>% count()
intersect(hart_locs, pcgr_locs) %>% count()
setdiff(hmf_locs, pcgr_locs) %>% count()
setdiff(pcgr_locs, hmf_locs) %>% count()















