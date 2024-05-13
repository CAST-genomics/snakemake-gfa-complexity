import sys
import os

from array import array


'''
Usage python plinklinear_to_bed.py [input.assoc.sorted] [output.bed]

'''
def main():
    inp, out = read_inputs()

    inp.readline()

    regions = get_sig_lines(inp)
    inp.close()

    filt_regions = filter_by_proximity(regions)

    print(f"Filtered out {len(regions)-len(filt_regions)} regions.")

    for line in filt_regions:
        out.write(line + "\n")

    out.close()

def read_inputs():
    # Verify inputs
    path_input = sys.argv[1]
    path_output = sys.argv[2]

    if not os.path.isfile(path_input):
        raise ValueError("Input .assoc path error!")

    inp = open(path_input, "r")
    out = open(path_output, "w")

    return inp, out

def get_sig_lines(inp, radius = 150_000):
    sig_lines = []
    for line in inp.readlines():
        # Beta corresponds to a z-score in this formulation.
        try:
            CHR, SNP, BP, A1, TEST, NMISS, BETA, STAT, P = line.split("\t")
            Z = float(BETA)
            if Z <= 5.326724:
                return sig_lines
            region = [CHR, max(int(BP) - radius, 0), int(BP) + radius, SNP, BETA]
            sig_lines.append(region)
        except Exception as e:
            print(line)
            raise e

def filter_by_proximity(regions, max_bins = 867, length_bin = 300_000):
    filtered_regions = []
    chrs = dict()
    for region in regions:
        CHR, START, END, SNP, SCORE = region
        start_bin = START // length_bin
        end_bin = END // length_bin
        cur_regions = chrs.setdefault(CHR, array('B', [0] * max_bins))
        if cur_regions[start_bin] or cur_regions[end_bin]:
            continue
        else:
            cur_regions[start_bin] = 1
            cur_regions[end_bin] = 1
            filtered_regions.append('\t'.join(str(item) for item in region))
    return filtered_regions


def snakemake_call(path_input, path_output):

    if not os.path.isfile(path_input):
        raise ValueError("Input .assoc path error!")

    inp = open(path_input, "r")
    out = open(path_output, "w")

    inp.readline()

    regions = get_sig_lines(inp)
    inp.close()

    filt_regions = filter_by_proximity(regions)

    print(f"Filtered out {len(regions)-len(filt_regions)} regions.")

    for line in filt_regions:
        out.write(line + "\n")

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

