
# Introduction

`rempipe` is a python library for running scientific workflows, or pipelines, such as RNA-Seq or SNP analysis. It has been designed to be modular, OS-independent, and most importantly, high-throughput.

> TODO: Add documentation instructions regarding usage and third-party install.

## Availability

`rempipe`, its source, and its documentation are licensed under the [GPL version 3](LICENSE.txt). The package does not contain any of third-party software required for running an analysis, e.g., tophat. All code and documentation included in `rempipe` was written by the authors, except where explicitly stated.

## Requirements

* Python v3.4.3 (or greater?)
* Desired third-party software, e.g., tophat.

# Installation

1. Open up the terminal ( <kbd>&#8984;</kbd> + <kbd>space</kbd> + `terminal` + <kbd>return</kbd> ) and change to the directory where you want to install.
![Terminal gif](https://dl.dropboxusercontent.com/u/1209050/screenshots/terminal.gif)

2. Clone the git repository and descend into the directory

    ```shell
git clone git@github.com:muppetjones/rempipe.git
cd rempipe
    ```

3. Run the unittests.

    ```shell
python -m unittest discover ./
    ```

# Available Pipelines

`rempipe` is designed to allow you to design and run your own pipelines; however, a number of preset pipelines exist for easier use.

>#### Initial QC :lock:
> - FastQC
> - Skewer
>
>Remove adapter sequences and trim reads using skewer. Read statistics are calculated before and after using FastQC

>#### Alignment
> - hisat2
> - samtools sort
> - samtools index
> - htseq-count
>
> Align reads to a reference genome, then count the number of reads aligned at each position.

>#### Assembly :lock:
> - velvet
> - abacas
> - RAST
>
> Assemble, sort, and annotate a bacterial genome. Includes several scripts for reporting assembly statistics.


>### RNA-Seq Pipeline
> - Initial QC pipeline
> - Alignment pipeline
> - DESeq2 analysis (via R) :lock:
>
> Handles alignment of reads to multiple genomes, including filtering reads through subsequent alignments. Also runs a differential expression analysis in `R` via a custom script.


>### WGS Pipeline :lock:
> - Initial QC pipeline
> - Assembly pipeline
> - SNP analysis (currently via kSNP3--kSNP4 coming soon)
>
> Everything you need from assembly and annotation to a k-mer based SNP analysis.

:lock: : Feature is currently unavailable on the `master` or `dev` branches.




# Testing and Known Bugs

`rempipe` has been tested on OSX 10.9, OSX 10.10, and CentOS with python 3.4.3.

* [Known Python 3.4-6 bug](https://bugs.python.org/issue18622)

    This bug is triggered in `tests.decorators.test_io`. The offending test
    is skipped by default, but may be reenabled after patching Python. The
    required patch is included in the `patch` directory.

    To apply the patch, run the following command:
    ```shell
    patch -p2 -b < file.patch
    ```
    where `file.patch` is the patch file.
