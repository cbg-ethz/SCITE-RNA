#!/bin/bash
#SBATCH --job-name=phylinsic
#SBATCH --output=logs/output_%A_%a.log
#SBATCH --error=logs/error_%A_%a.err
#SBATCH --array=0-27  # 28 entries, so 0-27
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=8G
#SBATCH --time=30:00:00

param_combinations=(
    "_param_testing/dropout_0"
    "_param_testing/dropout_0_2"
    "_param_testing/dropout_0_4"
    "_param_testing/dropout_0_6"
    "_param_testing/overdispersion_Het_3"
    "_param_testing/overdispersion_Het_6"
    "_param_testing/overdispersion_Het_10"
    "_param_testing/overdispersion_Het_100"
    "_param_testing/overdispersion_Hom_3"
    "_param_testing/overdispersion_Hom_6"
    "_param_testing/overdispersion_Hom_10"
    "_param_testing/overdispersion_Hom_100"
    "_param_testing/error_rate_0_001"
    "_param_testing/error_rate_0_01"
    "_param_testing/error_rate_0_05"
    "_param_testing/error_rate_0_1"
    "_param_testing/coverage_mean_10"
    "_param_testing/coverage_mean_30"
    "_param_testing/coverage_mean_60"
    "_param_testing/coverage_mean_100"
    "_param_testing/coverage_zero_inflation_0"
    "_param_testing/coverage_zero_inflation_0_2"
    "_param_testing/coverage_zero_inflation_0_4"
    "_param_testing/coverage_zero_inflation_0_6"
    "_param_testing/coverage_dispersion_1"
    "_param_testing/coverage_dispersion_2"
    "_param_testing/coverage_dispersion_5"
    "_param_testing/coverage_dispersion_10"
    "_param_testing/CNV_fraction_0_1"
    "_param_testing/CNV_fraction_0_3"
    "_param_testing/CNV_fraction_0_5"
    "_param_testing/homoplasy_fraction_0"
    "_param_testing/homoplasy_fraction_0_05"
    "_param_testing/homoplasy_fraction_0_1"
    "_param_testing/homoplasy_fraction_0_2"
)

APPENDIX=${param_combinations[$SLURM_ARRAY_TASK_ID]}

# Run Snakemake with the selected APPENDIX
snakemake --config APPENDIX="$APPENDIX" --use-conda -j 5 --rerun-incomplete
