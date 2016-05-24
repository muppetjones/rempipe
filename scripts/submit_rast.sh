#!/bin/bash -l

ROOT_DIR="$1"
cd "$ROOT_DIR"

# add RAST scripts to path
TOP="/Applications/myRAST.app";
if [[ "$PATH" != *${TOP}* ]]; then
    export PATH=$PATH:$TOP/bin
fi
if [[ "$PERL5LIB" != *${TOP}* ]]; then
    export PERL5LIB=$PERL5LIB:$TOP/lib:$TOP/modules/lib
fi

uname="muppetjones"
pswrd="israel43ongi"

# loop over all names in the genome names file
# NOTE: genome_names_sorted created from summarize_assembly script
#       -- should NOT include samples that were NOT assembled correctly
for name in $(cat "$ROOT_DIR"/genome_names_sorted.txt); do

    echo "Genome: $name"

    # loop over all of the contigs.fa files
    # -- essentially, any velvet assembly
    for assembly in $(find "$(pwd)" -type f -path "*${name}*" -name '*MULTI*.fa'); do

        echo "Assembly: $(basename $(dirname "$assembly"))/$(basename "$assembly")"

        odir="$(dirname "${assembly}")"
        svr_submit_RAST_job --user $uname --passwd $pswrd \
            --domain Bacteria --bioname "Mycobacterium tuberculosis" \
            --fasta "$assembly" 2>&1 | tee "${odir}"/rast.txt

        grep "Job" "${odir}"/rast.txt | sed -e "s/.*'\([0-9]\+\)'.*/\1/" > "${odir}"/rast.id
    done
done
