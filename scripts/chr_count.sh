#!/bin/bash

#' Count the number of reads that aligned to each chromosome.
#'
#'  Arguments:
#'      CHR_FILE    A file containing the string chromosome identifiers.
#'      BAM_FILE    BAM file(s) or a file pattern.
#'  Output:
#'      A count file for each BAM file, with extension '.ccnt'
#'      A count matrix in the *ASSUMED* common prefix directory
#'
#' NOTE: The script assumes that the BAM files are in sub directories.

USAGE="Usage: $0 <chr list file> <BAM FILE | pattern> [BAM FILE, ...]"
EXAMPLE="  EX $ $0 ~/genomes/human/human.chr ~/project/samples/*/*_human.bam"

# Read in the args--assume all but first are bam files!
# NOTE: Any files with spaces in the name MUST be single quoted
ARGS="$@"
CHR_FILE="$1"
shift
BAM_FILES="$@"

# Ensure the user gave a chromosome listing file
if [[ $CHR_FILE == *".bam" ]]; then
    echo "Chromosome list cannot be a BAM file"
    echo "$USAGE"
    echo "$EXAMPLE"
    exit 1
fi

# Parse the out matrix file name
FIRST_FILE="$(echo "$BAM_FILES" | cut -f 1 -d ' ')"
GENOME="$(basename ${CHR_FILE%%.*})"
MAT_OUT="$(dirname $(dirname "$FIRST_FILE"))/${GENOME}_chr.cov"
CHR_LIST="$(tail -n +2 "$CHR_FILE")"

# attempt to extract the number of base pairs (length) of each chromosome
SAM_HEADER="$(samtools view -H $FIRST_FILE | grep -E '^@.*SN:\w+\s' | sed -re  's/.*SN:(\w+)\s+LN:([0-9]+)/\1\t\2/')"
CHRLEN_FILE="${CHR_FILE}.len"
echo -e "chr\tlen\n$SAM_HEADER" > "$CHRLEN_FILE"

# inform the user to ensure expected data
echo "BAM files:"
for x in $BAM_FILES; do
    echo "   $(basename $x)";
done
echo -e "Chromosomes:\n   $(echo ${CHR_LIST} | sed -e 's/\s+/, /')"

# pause to allow user to cancel
# sleep 5

echo "Executing..."
for f in $BAM_FILES; do
    BAM_FILE=$f;
    OUT_FILE=${f/%bam/ccnt};

    if [ -f "$OUT_FILE" ]; then
        continue;
    fi
    echo "   $(basename $f)"

    # ensure the out file is a DIFFERENT name than the bam file
    # -- overwritting would be VERY, VERY bad
    if [[ "$BAM_FILE" == "$OUT_FILE" ]]; then
        echo 'Count will overwrite bam! $BAM_FILE > $OUTFILE'
        exit 1
    fi

    # put the sample name (assumed to be the last dir in path) as first line
    echo "$(basename $(dirname "$OUT_FILE"))" > $OUT_FILE;

    # count number of reads in each chromosome
    for n in $CHR_LIST; do
        samtools view $BAM_FILE $n | wc -l >> $OUT_FILE;
    done
done


# print the args used (for easier tracking)
echo "# $0 $ARGS" > "$MAT_OUT"
paste "$CHRLEN_FILE" \
    ${BAM_FILES//.bam/.ccnt} \
     >> "$MAT_OUT"

perl -pi -e "s/\t\s*/\t/g" ${MAT_OUT}
