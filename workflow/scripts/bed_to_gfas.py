import subprocess
import os
import sys

"""
For snakemake calls

Usage
python bed_to_gfas.py [bed_file] [out_dir]
"""
def main():

    bed_file, out_dir = read_inputs()
    gfab_file = "/home/wwford/datasets/hprc-v1.1-mc-grch38.gfab"

    gfabase_sub(gfab_file, bed_file, out_dir)

def read_inputs():
    bed_file = sys.argv[1]
    out_dir = sys.argv[2]

    check_inputs(bed_file, out_dir)

    return bed_file, out_dir

def check_inputs(bed_file, out_dir):
    if not os.path.isfile(bed_file):
        raise ValueError("Input bed file path error!")

    if not os.path.isdir(out_dir):
        print(f"{out_dir} does not exist, creating it")
        os.mkdir(out_dir)

def gfabase_sub(gfab_file, bed_file, out_dir):
    with open(bed_file, 'r') as bed:
        bed.readline()
        for segment in bed.readlines():
            vars = segment.strip().split("\t")
            print(vars)
            chr = vars[0]
            start = str(int(vars[1]) + 1)
            end = vars[2]

            out_gfa = os.path.join(out_dir, f"{chr}_{start}_{end}.gfa")

            command = ["gfabase","sub", gfab_file, "--range", 
                    f"GRCh38#{chr}:{start}-{end}", "--view", "-o", out_gfa]
            subprocess.run(command)

def snakemake_call(bed_file, out_file):

    out_dir = os.path.dirname(out_file)
    check_inputs(bed_file, out_dir)

    gfab_file = "/home/wwford/datasets/hprc-v1.1-mc-grch38.gfab"

    gfabase_sub(gfab_file, bed_file, out_dir)

    # Let snakemake know completion of all gfa generation
    open(os.path.join(out_dir,".success"), "w").close()

if __name__ == "__main__":
    try:
        import snakemake
        # Check if the script is executed within a Snakemake rule
        if snakemake:
            # Call my_function with input and output files from Snakemake
            snakemake_call(snakemake.input[0], snakemake.output[0])
    except ImportError:
        # Handle the case where snakemake module is not available
        main()

