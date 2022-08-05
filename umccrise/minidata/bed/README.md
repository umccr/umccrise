- `hg38_refseq_gencode_all_genes_v1.bed`:
  - "refseq hg38p.7 + gencodeV41 + recent refseq + 200 bp padding at both ends."
    (see [slack](https://umccr.slack.com/archives/C025TLC7D/p1659576393355259?thread_ts=1659521772.764189&cid=C025TLC7D))
  - for this copy, I kept only auto/sex contigs (removed those with `_` and `chrM`):
    - `cat hg38_refseq_gencode_all_genes_v1.bed | grep -E -v "_|chrM" | sort -V -k1,1 -k2,2`
    - `bgzip` then `tabix -p bed`

