#!/usr/bin/env bash

git pull
cd ngs_utils && git pull && cd ..
cd reference_data && git pull && cd ..
cd vcf_stuff && git pull && cd ..

