# Usage ./generate_node_sample_map.sh [dataset]

# Where $GFA refers to the minigraph cactus gfa for example
GFA=hprc-v1.1-mc-grch38.gfa.gz
FILE_PATH=hprc-v1.1-mc-grch38-node-sample-map-tabix.tsv.gz
zcat $GFA | grep -E '^W' | cut -f 2,7 | sed 's/\t//;s/</\t/g;s/>/\t/g' | awk -F $'\t' -v 'OFS=\t' '{
  for (i = 1; i <= NF; i++) {
    if (i == 1) {
      key = $i;
    } else {
      print($i, key);
    }
  }
}' | sort -u -t $'\t' -k1,1n -k2,2 | awk -F $'\t' -v 'OFS=\t' '
current==$1 { line = line OFS $2; next; }
{ print line; current=$1; line=$0; }
END { print line }' | awk -F $'\t' -v 'OFS=\t' 'NF {
  print "" FS $1 FS $0
  }'| bgzip > $FILE_PATH

# Compress file into tabix format for efficient indexing
tabix -s 1 -b 2 -e 3 $FILE_PATH

# Generate auxiliary file required for use of gfabase segmentation 
#   i.e. generating subgraphs
zcat hprc-v1.1-mc-grch38.gfa.gz | ./gfabase load -o hprc-v1.1-mc-grch38.gfab


# Pipe above output into below script to get tabix indexed node-haplotype map
#"""
#zcat $FILE_PATH | awk -F $'\t' -v 'OFS=\t' 'NF {print "" FS $1 FS $0}'| \
#       bgzip > minigraph_cactus_node_sample_map_tabix.tsv.gz
#
#tabix -s 1 -b 2 -e 3 minigraph_cactus_node_sample_map_tabix.tsv.gz

# Example query usage, note that segment name is an empty string """
#tabix minigraph_cactus_node_sample_map_tabix.tsv.gz :1-1

# Example File Output
#'''
#       1       1       HG00438 HG00621
#       2       2       HG00438
#       3       3       HG00438 HG01952
#       4       4       HG00438
#       5       5       HG00438
#       6       6       HG00438
#       7       7       HG00438
#       8       8       HG00621
#       9       9       HG00438 HG00621
#'''
#
#"""