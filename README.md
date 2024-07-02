# Snakemake workflow: GFA Complexity

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)

A Snakemake workflow for calculating the complexity of a GFA file.


## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).


## Notes:

The snakemake.old file was used to generate several regions, but it requires a good amount of manual work to prepare each run. So I tried to make the snakemake file but it hangs when running on slurm nodes on expanse (but not on the login node for dry-runs) and I'm not sure why.

Our complexity metric scales with the size of the genome so I had been using a random background sequence to compare against and determine signifcance. Though it requires more testing, I think using a highly polygenic trait's `sumstats` file like height gwas hits actually produces a more realistic background distribution.

Only minigraph cactus version 1.1 of the pangenome has the walk labels we use to generate the node-hap map so we are limited to this version. It uses GRCh38 chromosomes and labels to build off of.

## Workflow

All scripts should be in `workflow/scripts/`

1. Generate a node-haplotype map using `build_node_sample_map.sh` and build a gfabase `.gfab` file for fast graph subsetting. Note: this file is not an executable and more a collection of one liners that I saved in a single location.
2. Build `bed` files for graph segmenting. From `.sumstats` file formats I built the scripts `sumstats_to_plinklinear.py` and `plinklinear_to_bed.py` to generate bed files. It requires a snp to GRCh38 coordinate file I currated at `/home/wwford/datasets/rsID.coords.bed` on Expanse. Anyone should be able to read this file.
3. Segment subgraphs using `bed_to_gfas.py`.
4. Calculate complexity using `gfa_to_complexity.py`

There is another script `manhattan_range.py` adapted from Arya's script that takes in plinklinear files and plots manhattan plots to verify that we are accurately grabbing the significant regions.

## Software Requirements

In the workflow above we use that will likely need to be downloaded if not already:
GFA Base
tabix
bgzip

snakemake (though I didn't get the efficient version of this to work)


rsID.coords.bed file example:
```
chr10	132889560	132889561	rs11815293
chr2	156133793	156133794	rs12694551
...
```

Feel free to email me if you have any questions: wwford@ucsd.edu

