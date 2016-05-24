#!/bin/bash
# USAGE:
#   ./summarize_assembly.sh <dir> [pattern]
#
# Example:
#   ./summarize_assembly.sh ~/data/project/samples '/*/assembly/contigs.fa'
#
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
if [ -z ${2+x} ]; then
    pattern="/*/assembly/contigs.fa"  # velvet contigs
    pattern="/*/assembly/*MULTIFASTA.fa"  # abacas contigs
else
    pattern="$2"
fi

echo -e "Using...\n\tdirectory '$absin'\n\tpattern '$pattern'"

# genome list
genomelist="${absin}/genome_list_all.txt"
rm -f "$genomelist"
for f in "${absin}"${pattern}; do
    # skip if the file is empty (-s)
    if [ -s "$f" ]; then
        genomedir="$(dirname "$(dirname "$f")")"
        genomename="$(basename "${genomedir}")"
        echo -e "${f}\t${genomename}" >> "${genomelist}"
    fi
done

### annotated geneome list
### -- this will be a list of the genome names, sorted by the longest contigs

# create file of contig sizes for each genome
for f in ${absin}${pattern}; do
    # skip if the file is empty (! -s)
    if [ -s "$f" ]; then
        contigsize="$(dirname "$f")/contig_sizes.txt"
        grep ">" "$f" |awk -F '_' '{print $4}' |sort -rn > "${contigsize}"
    fi
done

# genome list
genomelist="${absin}/genome_list.txt"
rm -f "$genomelist"
for f in "${absin}"${pattern}; do
    # skip if the file is empty (-s)
    if [ -s "$(dirname "$f")/contig_sizes.txt" ]; then
        genomedir="$(dirname "$(dirname "$f")")"
        genomename="$(basename "${genomedir}")"
        echo -e "${f}\t${genomename}" >> "${genomelist}"
    fi
done

# create a list of empty contig size files
empty="${absin}/genome_list_empty.txt"
echo "# $(date)" > "${empty}"
for f in "${absin}"/*/assembly/contig_sizes.txt; do
    # skip if the file is NOT empty (! -s)
    if [ ! -s "$f" ]; then
        genomename="$(basename "$(dirname "$(dirname "$f")")")"
        echo "${genomename}" >> "${empty}"
    fi
done

# create a summary file of all assemblies
summary="${absin}/assembly_summary.txt"
# rm -f "$summary"
echo "# $(date)" > "${summary}"
echo "# $0 $@" >> "${summary}"
for f in ${absin}/*/assembly/contig_sizes.txt; do
    echo "$f"
    genomename="$(basename "$(dirname "$(dirname "$f")")")"
    if [ -s "$f" ]; then
        num_contig="$(wc -l < "$f")"
        max_contig_size="$(head -n1 "$f")"
        total_bp="$(perl -nle '$sum += $_ } END { print $sum' "$f")"
    else
        num_contig='       0'  # strange spacing due to cmd?
        max_contig_size='     -'
        total_bp='     -'
    fi
    echo -e "${genomename}\t${num_contig}\t${max_contig_size}\t${total_bp}" \
        >> "$summary"
done

# create a sorted list (by contig size, largest first) of the genome names
# -- remove any comment lines and any lines with '-'
sortednames="${absin}/genome_names_sorted.txt"
cat "$summary" | grep -v "^#" | awk '$3 == "-" {next} {print}' \
    | sort -k3 -rn | awk '{print $1}' > "$sortednames"
