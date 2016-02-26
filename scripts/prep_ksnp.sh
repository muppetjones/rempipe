#!/bin/bash -l
# This script summarizes a set of genome assemblies and calculates the
# optimal k to use for kSNP3.
# NOTE: Do NOT NOT NOT run this on qsub. Kchooser calls jellyfish via
#       a system call that makes an inappropriate ioctl call.

set -e

# scripts and programs
ksnp_dir="${HOME}/dev/kSNP3_Mac/kSNP3"

# get absolute directory paths (and make sure odir exists)
# idir="$(readlink -f "${1}")"
idir="${1}"
cd "${idir}"

echo "$idir"

# set filenames
genome_list="${idir}/genome_list.txt"
genome_names="${idir}/genome_names_sorted.txt"
kreport="${idir}/Kchooser.report"

# run 'summarize_assembly' if necessary
if [ ! -f "$genome_list" ]; then
    echo "-----------------------"
    echo "Summarizing assemblies in '${idir}'..."
    ~/dev/pipe/scripts/summarize_assembly.sh "${idir}"
    echo "-----------------------"
    echo "Created files:"
    echo "  genome list:    ${genome_list}"
    echo "  genome names:   ${genome_names}"
fi

echo $PATH

# run Kchooser if necessary
# (kind of...run the script that sets up and runs Kchooser)
echo "$ksnp_dir"
echo "----------------------"
echo "Calculating optimal k for assemblies in '${idir}'"
set -x
"${ksnp_dir}"/MakeFasta "${genome_list}" "${idir}all_contigs.fa" <<EOL
C
EOL
"${ksnp_dir}"/Kchooser "${idir}all_contigs.fa"
set +x
rm -f "${idir}/all_contigs.fa"
echo "----------------------"
k="$(grep "optimum" Kchooser.report | cut -d' ' -f7 | cut -d'.' -f1)"
ncpu=8

# Quick summary
echo "Parameters:"
echo "  k:              ${k}"
echo "  ncpu:           ${ncpu}"
echo "------------------------------------------------------"
