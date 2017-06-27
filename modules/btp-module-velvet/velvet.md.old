Key Learning Outcomes
---------------------

After completing this practical the trainee should be able to:

-   Compile velvet with appropriate compile-time parameters set for a
    specific analysis

-   Be able to choose appropriate assembly parameters

-   Assemble a set of single-ended reads

-   Assemble a set of paired-end reads from a single insert-size library

-   Be able to visualise an assembly in AMOS Hawkeye

-   Understand the importance of using paired-end libraries in *de novo*
    genome assembly

Resources You’ll be Using
-------------------------

Although we have provided you with an environment which contains all the
tools and data you will be using in this module, you may like to know
where we have sourced those tools and data from.

### Tools Used

Velvet:  
http://www.ebi.ac.uk/~zerbino/velvet/

AMOS Hawkeye:  
http://apps.sourceforge.net/mediawiki/amos/index.php?title=Hawkeye

gnx-tools:  
https://github.com/mh11/gnx-tools

FastQC:  
http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/

R:  
http://www.r-project.org/

### Sources of Data

-  ftp://ftp.ensemblgenomes.org/pub/release-8/bacteria/fasta/Staphylococcus/s_aureus_mrsa252/dna/s_aureus_mrsa252.EB1_s_aureus_mrsa252.dna.chromosome.Chromosome.fa.gz

-   http://www.ebi.ac.uk/ena/data/view/SRS004748

-   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR022/SRR022825/SRR022825.fastq.gz

-   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR022/SRR022823/SRR022823.fastq.gz

-   http://www.ebi.ac.uk/ena/data/view/SRX008042

-   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR022/SRR022852/SRR022852_1.fastq.gz

-   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR022/SRR022852/SRR022852_2.fastq.gz

-   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR023/SRR023408/SRR023408_1.fastq.gz

-   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR023/SRR023408/SRR023408_2.fastq.gz

-   http://www.ebi.ac.uk/ena/data/view/SRX000181

-   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR000/SRR000892/SRR000892.fastq.gz

-   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR000/SRR000893/SRR000893.fastq.gz

-   http://www.ebi.ac.uk/ena/data/view/SRX007709

-   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR022/SRR022863/SRR022863_1.fastq.gz

-   ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR022/SRR022863/SRR022863_2.fastq.gz

Introduction
------------

The aim of this module is to become familiar with performing *de novo*
genome assembly using Velvet, a de Bruijn graph based assembler, on a
variety of sequence data.

Prepare the Environment
-----------------------

The first exercise should get you a little more comfortable with the
computer environment and the command line.

First make sure that you are in the denovo working directory by typing:

    cd /home/trainee/denovo

and making absolutely sure you’re there by typing:

    pwd

Now create sub-directories for this and the two other velvet practicals.
All these directories will be made as sub-directories of a directory for
the whole course called NGS. For this you can use the following
commands:

    mkdir -p NGS/velvet/{part1,part2,part3}

The `-p` tells `mkdir` (make directory) to make any parent directories
if they don’t already exist. You could have created the above
directories one-at-a-time by doing this instead:

    mkdir NGS
    mkdir NGS/velvet
    mkdir NGS/velvet/part1
    mkdir NGS/velvet/part2
    mkdir NGS/velvet/part3

After creating the directories, examine the structure and move into the
directory ready for the first velvet exercise by typing:

    ls -R NGS
    cd NGS/velvet/part1
    pwd

Downloading and Compiling Velvet
--------------------------------

For the duration of this workshop, all the software you require has been
set up for you already. This might not be the case when you return to
“real life”. Many of the programs you will need, including velvet, are
quite easy to set up, it might be instructive to try a couple.

Although you will be using the preinstalled version of velvet, it is
useful to know how to compile velvet as some of the parameters you might
like to control can only be set at compile time. You can find the latest
version of velvet at:

http://www.ebi.ac.uk/~zerbino/velvet/

You could go to this URL and download the latest velvet version, or
equivalently, you could type the following, which will download, unpack,
inspect, compile and execute your locally compiled version of velvet:

    cd /home/trainee/denovo/NGS/velvet/part1
    pwd
    tar xzf /home/trainee/denovo/data/velvet_1.2.10.tgz
    ls -R
    cd velvet_1.2.10
    make
    ./velveth

The standout displayed to screen when ’make’ runs may contain an error
message but it is ignored

Take a look at the executables you have created. They will be displayed
as green by the command:

    ls --color=always

The switch `–color`, instructs that files be coloured according to their
type. This is often the default but we are just being explicit. By
specifying the value `always`, we ensure that colouring is always
applied, even from a script.

Have a look of the output the command produces and you will see that
`MAXKMERLENGTH=31` and `CATEGORIES=2` parameters were passed into the
compiler.

This indicates that the default compilation was set for de Bruijn graph
k-mers of maximum size 31 and to allow a maximum of just 2 read
categories. You can override these, and other, default configuration
choices using command line parameters. Assume, you want to run velvet
with a k-mer length of 41 using 3 categories, velvet needs to be
recompiled to enable this functionality by typing:

    make clean
    make MAXKMERLENGTH=41 CATEGORIES=3
    ./velveth

Discuss with the persons next to you the following questions:\
What are the consequences of the parameters you have given make for
velvet?

MAXKMERLENGTH: increase the max k-mer length from 31 to 41

CATEGORIES: paired-end data require to be put into separate categories.
By increasing this parameter from 2 to 3 allows you to process 3 paired
/ mate-pair libraries and unpaired data.

Why does Velvet use k-mer 31 and 2 categories as default?

Possibly a number of reason:\
- odd number to avoid palindromes\
- The first reads were very short (20-40 bp) and there were hardly any
paired-end data\
around so there was no need to allow for longer k-mer lengths / more
categories.\
- For programmers: 31 bp get stored in 64 bits (using 2bit encoding)

Should you get better results by using a longer k-mer length?

If you can achieve a good k-mer coverage - yes.

What effect would the following compile-time parameters have on velvet:\
`OPENMP=Y`

Turn on multithreading

`LONGSEQUENCES=Y`

Assembling reads / contigs longer than 32kb long

`BIGASSEMBLY=Y`

Using more than 2.2 billion reads

`VBIGASSEMBLY=Y`

Not documented yet

`SINGLE_COV_CAT=Y`

Merge all coverage statistics into a single variable - save memory

For a further description of velvet compile and runtime parameters
please see the velvet Manual:
https://github.com/dzerbino/velvet/wiki/Manual
