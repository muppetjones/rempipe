#!/bin/bash
# Summarize velvet assemblies.
#
# Creates files intended for use in kSNP3, but just useful for summarizing
# a set of velvet assemblies.
#   1) genome_list.txt -- a list of contig paths with the genome name (tab)
#   2) */contig_sizes.txt -- a list of sorted contig sizes for each genome
#   3) assembly_summary.txt -- genome_name  n_contigs  length_largest_contig
#   4) assembly_names_sorted.txt -- genome names sorted by contig size
#
# Examples:
#   1) genome_list.txt -- unordered (should be by name)
#       /path/to/data/assembly/genome1/contigs.fa   genome1
#           ...
#       /path/to/data/assembly/genomeN/contigs.fa   genomeN
#   2) genome1/contig_sizes.txt -- one for every genome!
#       300300
#       280100
#       ...
#       412
#   3) assembly_summary.txt -- unordered (should be by name)
#       genome1     119     300300
#       ...
#       genomeN     1395    4610
#   4) assembly_names_sorted.txt -- ordered largest to smallest contig size
#       genome6
#       genome1
#       ...
#       genomeN
#
# NOTE: As stated in the kSNP3 guide--you cannot move these files!
#
# NOTE: It is assumed that the given directory contains all of the directories
#       given to velvetg, i.e., all genomes are found in ./*/contigs.fa

yell() { echo -e "\n$0: $*" >&2; }
die() { yell "ERROR: $*"; exit 111; }
try() { "$@" || die "cannot $*"; }

if [ -z $1 ]; then
    die "No input directory given"
fi

indir="$1"
absin="$(cd "$indir"; pwd)"

# genome list
genomelist="${absin}/genome_list.txt"
rm -f "$genomelist"
for f in "${absin}"/*/assembly/contigs.fa; do
    genomedir="$(dirname "$(dirname "$f")")"
    genomename="$(basename "${genomedir}")"
    echo -e "${f}\t${genomename}" >> "${genomelist}"
done

### annotated geneome list
### -- this will be a list of the genome names, sorted by the longest contigs

# create file of contig sizes for each genome
for f in "${absin}"/*/assembly/contigs.fa; do
    contigsize="$(dirname "$f")/contig_sizes.txt"
    grep ">" "$f" |awk -F '_' '{print $4}' |sort -rn > "${contigsize}"
done

# create a summary file of all assemblies
summary="${absin}/assembly_summary.txt"
rm -f "$summary"
for f in "${absin}"/*/assembly/contig_sizes.txt; do
    genomename="$(basename "$(dirname "$(dirname "$f")")")"
    echo -e "${genomename}\t$(wc -l < "$f")\t$(head -n1 "$f")" >> "$summary"
done

# create a sorted list (by contig size, largest first) of the genome names
sortednames="${absin}/genome_names_sorted.txt"
sort -k3 -rn "$summary" | awk '{print $1}' > "${sortednames}"
