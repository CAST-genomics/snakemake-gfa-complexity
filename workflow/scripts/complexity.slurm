#!/usr/bin/bash
#SBATCH --export ALL
#SBATCH --partition ind-shared
#SBATCH --account ddp268
#SBATCH --job-name test
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --time 00:10:00
#SBATCH --output /home/wwford/analysis/UKB_460K.body_HEIGHTz/slurm/%J.out-%N
#SBATCH --error /home/wwford/analysis/UKB_460K.body_HEIGHTz/slurm/%J.err-%N

#e.g. sbatch -J test complexity.slurm chr1_874462_1174461.gfa test.score

WD="/home/wwford/analysis/UKB_460K.body_HEIGHTz/snakemake-gfa-complexity/workflow/scripts/"
INDIR="/home/wwford/analysis/UKB_460K.body_HEIGHTz/subset_graphs_10e-8/"
OUTDIR="/home/wwford/analysis/UKB_460K.body_HEIGHTz/scores/"

echo "$INDIR""$1" "$OUTDIR""$2"
python "$WD"gfa_to_complexity.py "$INDIR""$1" "$OUTDIR""$2"
