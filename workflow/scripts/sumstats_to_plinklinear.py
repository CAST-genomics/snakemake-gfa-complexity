#.sumstats file format ex:
#```
#SNP     A1      A2      N       Z
#rs3094315       A       G       14916.610       -0.253
#rs1048488       T       C       14481.670       -0.238
#```
#Corresponding .assoc.linear file format ex:
#```
#CHR SNP BP  A1  TEST    NMISS   BETA   STAT    P
#M   rs3094315   00  A   file_name   14916.610   *  -0.253  000
#M   rs1048488   00  T   file_name   14481.670   *  -0.238  000
#```

# Usage: python sumstats_to_plinklinear.py <trait.sumstats> <output>
# i.e. python sumstats_to_plinklinear.py /expanse/lustre/projects/ddp412/tamariutabartell/sumstats/UKB_460K.body_HEIGHTz.sumstats UKB_460K.body_HEIGHTz.assoc

import subprocess
import sys
import os

import pandas as pd
from scipy.stats import norm

def main():
    inp, out = read_inputs()

    header = "\t".join([
        "#CHROM", "ID", "POS", "A1", "TEST",
        "NMISS", "BETA", "STAT", "P"
    ])+"\n"
    out.write(header)

    mapper = coordsMapper()
    inp.readline()
    for line in inp.readlines():
        formatted_line = format_line(line, mapper)
        out.write(formatted_line+"\n")

    inp.close()
    out.close()


def read_inputs():
    # Verify inputs
    path_input = sys.argv[1]
    path_output = sys.argv[2]

    if not os.path.isfile(path_input):
        raise ValueError("Input sumstats path error!")

    inp = open(path_input, "r")
    out = open(path_output, "w")

    return inp, out

def format_line(line, mapper):
    line = line.strip()
    values = line.split("\t")
    if len(values) == 6:
        SNP, A1, A2, N, CHISQ, Z = values
    if len(values) == 5:
        SNP, A1, A2, N, Z = values
        CHISQ = "*"
    else:
        print(values)
        raise ValueError("Wrong number of input columns in sumstats file!")

    CHR, BP = mapper.get_coords(SNP)

    TEST = "*"
    NMISS = N
    STAT = CHISQ
    BETA = abs(float(Z))
    P = get_p(float(Z))

    values = [str(CHR), str(SNP), str(BP), str(A1), str(TEST), str(NMISS), str(BETA), str(STAT), str(P)]

    return "\t".join(values)

class coordsMapper():
    def __init__(self):
        coords_file = "~/datasets/rsID.coords.bed"
        coords_df = pd.read_csv(coords_file, sep="\t", 
            names = ["chr","pos","end","rsid"]).dropna()

        self.rsID_chrom = pd.Series(coords_df.chr.values, index=coords_df.rsid).to_dict()
        self.rsID_pos = pd.Series(coords_df.pos.values, index=coords_df.rsid).to_dict()
        self.leftovers = open("/home/wwford/datasets/leftover_rsIDs", "w+")

    def get_coords(self, SNP):
        try:
            return self.rsID_chrom[SNP], self.rsID_pos[SNP]
        except:
            self.leftovers.write(SNP)
            return 0, 0
    
    def __del__(self):
        self.leftovers.close()

def get_p(Z):
    return 1 - norm.cdf(abs(Z))

def snakemake_call(path_input, path_output):
    # Verify inputs

    if not os.path.isfile(path_input):
        raise ValueError("Input sumstats path error!")

    inp = open(path_input, "r")

    if not os.path.isdir(os.path.dirname(path_output)):
        os.mkdir(os.path.dirname(path_output))
    
    out = open(path_output, "w")

    header = "\t".join([
        "CHR", "SNP", "BP", "A1", "TEST",
        "NMISS", "BETA", "STAT", "P"
    ])+"\n"
    out.write(header)

    mapper = coordsMapper()
    inp.readline()
    for line in inp.readlines():
        formatted_line = format_line(line, mapper)
        out.write(formatted_line+"\n")

    inp.close()
    out.close()


if __name__ == "__main__":
    try:
        import snakemake
        # Check if the script is executed within a Snakemake rule
        if snakemake:
            # Call my_function with input and output files from Snakemake
            snakemake_call(snakemake.input[0], snakemake.output[0])
            return
    except ImportError:
        # Handle the non snakemake call
        main()
