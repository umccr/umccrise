#!/usr/bin/env bash

VCF=MB_100vs50-vardict-annotated.vcf.gz
#VCF=test.vcf.gz

OUT=${VCF/.vcf.gz/.cyvcf2_writer.vcf}
time zsh -c "python vardict/test_parsers/test_vcf_parsers.py ${VCF} cyvcf2 ${OUT} && bgzip -f ${OUT}"
# 291.77s user 6.91s system 92%  cpu 5:23.91 total

time zsh -c "python vardict/test_parsers/test_vcf_parsers.py ${VCF} cyvcf2 | bgzip -c > ${VCF/.vcf/.cyvcf2.vcf}"
# 283.23s user 2.40s system 183% cpu 2:35.55 total

time zsh -c "python vardict/test_parsers/test_vcf_parsers.py ${VCF} pysam | bgzip -c > ${VCF/.vcf/.pysam.vcf}"
# 325.09s user 2.67s system 171% cpu 3:11.47 total

time zsh -c "python vardict/test_parsers/test_vcf_parsers.py ${VCF} python | bgzip -c > ${VCF/.vcf/.python.vcf}"
# 251.17s user 1.84s system 179% cpu 2:20.60 total

time zsh -c "gunzip -c ${VCF} | py -x 'test_vcf_parsers.proc_line(x)' | bgzip -c > ${VCF/.vcf/.py.vcf}"
# 309.34s user 4.74s system 186% cpu 2:48.33 total




