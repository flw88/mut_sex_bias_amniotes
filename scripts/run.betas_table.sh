#!/usr/bin/env bash

sbatch --account="palab" -t 01:00:00 -c 1 --mem-per-cpu=12G -J betas_table -o out/betas_table ./betas_table.R
