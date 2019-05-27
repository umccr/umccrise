#!/usr/bin/env bash

git pull
cd ngs_utils && git pull && cd ..
cd hpc_utils && git pull && cd ..
cd vcf_stuff && git pull && cd ..

