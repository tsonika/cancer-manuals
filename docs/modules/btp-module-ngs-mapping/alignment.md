Key Learning Outcomes
---------------------

After completing this practical the trainee should be able to:

-   Perform a simple NGS data alignment task, with Bowtie2, against one
    interested reference data

-   Interpret and manipulate the mapping output using SAMtools

-   Visualise the alignment via a standard genome browser, e.g. IGV
    browser

Resources You’ll be Using
-------------------------

### Tools Used

Bowtie2:  
http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

Samtools:  
http://broadinstitute.github.io/picard

BEDTools:  
http://code.google.com/p/bedtools/

UCSC tools:  
http://hgdownload.cse.ucsc.edu/admin/exe/

IGV genome browser:  
http://www.broadinstitute.org/igv/

Useful Links
------------

SAM Specification:  
http://samtools.sourceforge.net/SAM1.pdf

Explain SAM Flags:  
https://broadinstitute.github.io/picard/explain-flags.html

### Sources of Data

http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-11431

## Author Information

*Primary Author(s):*
    Myrto Kostadima kostadim@ebi.ac.uk
   
*Contributor(s):*
    Xi (Sean) Li sean.li@csiro.au

Introduction
------------

The goal of this hands-on session is to perform an unspliced alignment
for a small subset of raw reads. We will align raw sequencing data to
the mouse genome using Bowtie2 and then we will manipulate the SAM
output in order to visualize the alignment on the IGV browser.

Prepare the Environment
-----------------------

We will use one data set in this practical, which can be found in the
`ChIP-seq` directory on your desktop.

Open the Terminal.

First, go to the right folder, where the data are stored.

    cd /home/trainee/chipseq

The `.fastq` file that we will align is called `Oct4.fastq`. This file
is based on Oct4 ChIP-seq data published by Chen *et al.* (2008). For
the sake of time, we will align these reads to a single mouse
chromosome.

Alignment
---------

You already know that there are a number of competing tools for short
read alignment, each with its own set of strengths, weaknesses, and
caveats. Here we will try Bowtie2, a widely used ultrafast, memory
efficient short read aligner.

Bowtie2 has a number of parameters in order to perform the alignment. To
view them all type

    bowtie2 --help

Bowtie2 uses indexed genome for the alignment in order to keep its
memory footprint small. Because of time constraints we will build the
index only for one chromosome of the mouse genome. For this we need the
chromosome sequence in FASTA format. This is stored in a file named
`mm10`, under the subdirectory `bowtie_index`.

The indexed chromosome is generated using the command:
!!! failure "STOP"
    DO NOT run this command. This has already been run for you.

    ** bowtie2-build bowtie_index/mm10.fa bowtie_index/mm10 **

This command will output 6 files that constitute the index. These files
that have the prefix `mm10` are stored in the `bowtie_index`
subdirectory. To view if they files have been successfully created type:

    ls -l bowtie_index

Now that the genome is indexed we can move on to the actual alignment.
The first argument for `bowtie2` is the basename of the index for the
genome to be searched; in our case this is `mm10`. We also want to make
sure that the output is in SAM format using the `-S` parameter. The last
argument is the name of the FASTQ file.

Align the Oct4 reads using Bowtie2:

    bowtie2 -x bowtie_index/mm10 -q Oct4.fastq > Oct4.sam

The above command outputs the alignment in SAM format and stores them in
the file `Oct4.sam`.

In general before you run Bowtie2, you have to know what quality
encoding your FASTQ files are in. The available FASTQ encodings for
bowtie are:

–phred33-quals
:   Input qualities are Phred+33 (default).

–phred64-quals
:   Input qualities are Phred+64 (same as `–solexa1.3-quals`).

–solexa-quals
:   Input qualities are from GA Pipeline ver. < 1.3.

–solexa1.3-quals
:   Input qualities are from GA Pipeline ver. >= 1.3.

–integer-quals
:   Qualities are given as space-separated integers (not ASCII).

    
The FASTQ files we are working with are Sanger encoded (Phred+33), which
is the default for Bowtie2.

Bowtie2 will take 2-3 minutes to align the file. This is fast compared
to other aligners which sacrifice some speed to obtain higher
sensitivity.

Look at the top 10 lines of the SAM file using head (record lines are
wrapped). Then try the second command, note use arrow navigation and to
exit type ’q’.

    head  Oct4.sam
    less -S Oct4.sam

!!! note "Question"
    Can you distinguish between the header of the SAM format and the actual alignments?

??? "**Answer**"
    !!! success "Answer"
        The header line starts with the letter ‘@’, i.e.:

        ----- ------------ -------------- ---------- ----------------------------------------------------------------------------------------------------------
        @HD   VN:1.0       SO:unsorted               
        @SQ   SN:chr1      LN:195471971              
        @PG   ID:Bowtie2   PN:bowtie2     VN:2.2.4   CL:“/tools/bowtie2/bowtie2-default/bowtie2-align-s –wrapper basic-0 -x bowtie\_index/mm10 -q Oct4.fastq”
        ----- ------------ -------------- ---------- ----------------------------------------------------------------------------------------------------------

While, the actual alignments start with read id, i.e.:

  -------------- ---- ------ -----
  SRR002012.45   0    etc    
  SRR002012.48   16   chr1   etc
  -------------- ---- ------ -----

!!! note "Question"
    What kind of information does the header provide?

??? "**Answer**"
    !!! success "Answer"
    
    -   @HD: Header line; VN: Format version; SO: the sort order of alignments.

    -   @SQ: Reference sequence information; SN: reference sequence name; LN: reference sequence length.

     -   @PG: Read group information; ID: Read group identifier; VN: Program version; CL: the command line that produces the alignment.

!!! note "Question"
    To which chromosome are the reads mapped?

??? "**Answer**"
    !!! success "Answer"
        Chromosome 1.

Manipulate SAM output
---------------------

SAM files are rather big and when dealing with a high volume of NGS
data, storage space can become an issue. As we have already seen, we can
convert SAM to BAM files (their binary equivalent that are not human
readable) that occupy much less space.

Convert SAM to BAM using `samtools view` and store the output in the
file `Oct4.bam`. You have to instruct `samtools view` that the input is
in SAM format (`-S`), the output should be in BAM format (`-b`) and that
you want the output to be stored in the file specified by the `-o`
option:

    samtools view -bSo Oct4.bam Oct4.sam

Compute summary stats for the Flag values associated with the alignments
using:

    samtools flagstat Oct4.bam

Visualize alignments in IGV
---------------------------

IGV is a stand-alone genome browser. Please check their website
(<http://www.broadinstitute.org/igv/>) for all the formats that IGV can
display. For our visualization purposes we will use the BAM and bigWig
formats.

When uploading a BAM file into the genome browser, the browser will look
for the index of the BAM file in the same folder where the BAM files is.
The index file should have the same name as the BAM file and the suffix
`.bai`. Finally, to create the index of a BAM file you need to make sure
that the file is sorted according to chromosomal coordinates.

Sort alignments according to chromosomal position and store the result
in the file with the prefix `Oct4.sorted`:

    samtools sort Oct4.bam -o Oct4.sorted.bam

Index the sorted file.

    samtools index Oct4.sorted.bam

The indexing will create a file called `Oct4.sorted.bam.bai`. Note that
you don’t have to specify the name of the index file when running
`samtools index`, it simply appends a `.bai` suffix to the input BAM
file.

Another way to visualize the alignments is to convert the BAM file into
a bigWig file. The bigWig format is for display of dense, continuous
data and the data will be displayed as a graph. The resulting bigWig
files are in an indexed binary format.

The BAM to bigWig conversion takes place in two steps. Firstly, we
convert the BAM file into a bedgraph, called `Oct4.bedgraph`, using the
tool `genomeCoverageBed` from BEDTools. Then we convert the bedgraph
into a bigWig binary file called `Oct4.bw`, using `bedGraphToBigWig`
from the UCSC tools:

    genomeCoverageBed -bg -ibam Oct4.sorted.bam -g bowtie_index/mouse.mm10.genome > Oct4.bedgraph
    bedGraphToBigWig Oct4.bedgraph bowtie_index/mouse.mm10.genome Oct4.bw

Both of the commands above take as input a file called
`mouse.mm10.genome` that is stored under the subdirectory
`bowtie_index`. These genome files are tab-delimited and describe the
size of the chromosomes for the organism of interest. When using the
UCSC Genome Browser, Ensembl, or Galaxy, you typically indicate which
species/genome build you are working with. The way you do this for
BEDTools is to create a “genome” file, which simply lists the names of
the chromosomes (or scaffolds, etc.) and their size (in basepairs).

BEDTools includes pre-defined genome files for human and mouse in the
`genomes` subdirectory included in the BEDTools distribution.

Now we will load the data into the IGV browser for visualization. In
order to launch IGV double click on the `IGV 2.3` icon on your Desktop.
Ignore any warnings and when it opens you have to load the genome of
interest.

On the top left of your screen choose from the drop down menu
`Mouse (mm10)`. If it doesn’t appear in list, click `More ..`, type
`mm10` in the Filter section, choose the mouse genome and press OK. Then
in order to load the desire files go to:

    File > Load from File

On the pop up window navigate to Desktop -> chipseq folder and select
the file `Oct4.sorted.bam`.

Repeat these steps in order to load `Oct4.bw` as well.

Select `chr1` from the drop down menu on the top left. Right click on
the name of `Oct4.bw` and choose Maximum under the Windowing Function.
Right click again and select Autoscale.

In order to see the aligned reads of the BAM file, you need to zoom in
to a specific region. For example, look for gene `Lemd1` in the search
box.

!!! note "Question"
    What is the main difference between the visualization of BAM and bigWig files?

??? "**Answer**"
    !!! success "Answer"
        The actual alignment of reads that stack to a particular region can be displayed using the information stored in a BAM format. The bigWig format is for display of dense, continuous data that will be displayed in the Genome Browser as a graph.

Using the `+` button on the top right, zoom in to see more of the
details of the alignments.

!!! note "Question"
    What do you think the different colors mean?

??? "**Answer**"
    !!! success "Answer"
        The different color represents four nucleotides, e.g. blue is Cytidine (C), red is Thymidine (T).

Practice Makes Perfect!
-----------------------

In the chipseq folder you will find the file `gfp.fastq`. Follow the
above described analysis, from the bowtie2 alignment step, for this
dataset as well. You will need these files for the ChIP-Seq module.
