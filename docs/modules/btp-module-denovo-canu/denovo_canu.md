# Pacbio reads: assembly with command line tools

<!-- *Keywords: de novo assembly, PacBio, PacificBiosciences, Illumina, command line, Canu, Circlator, BWA, Spades, Pilon, Microbial Genomics Virtual Laboratory* -->

This tutorial demonstrates how to use long PacBio sequence reads to assemble a bacterial genome, including correcting the assembly with short Illumina reads.

## Resources

Tools (and versions) used in this tutorial include:

- canu 1.5 (requires java 1.8)
- infoseq and sizeseq (part of EMBOSS) 6.6.0.0
- circlator 1.5.1
- bwa 0.7.15
- samtools 1.3.1
- makeblastdb and blastn (part of blast) 2.4.0+
- pilon 1.20
- spades 3.10.1
<!-- - prokka 1.12 -->

## Learning objectives

At the end of this tutorial, be able to:

1. Assemble and circularise a bacterial genome from PacBio sequence data.
2. Recover small plasmids missed by long read sequencing, using Illumina data
3. Explore the effect of polishing assembled sequences with a different data set.

## Overview

Simplified version of workflow:

![workflow](images/flowchart.png)

!!! warning "Important" 
    Background Information 
    
!!! note "Question"
        How do long- and short-read assembly methods differ?

!!! success ""
    ??? "**Answer**"
        Short reads are usually assembled using De Bruijn graphs. With long reads, there is a move back towards simpler overlap-layout-consensus methods.

!!! note "Question"
    Where can we find out the what the approximate genome size should be for the species being assembled?

!!! success ""
    ??? "**Answer**"
        Go to https://www.ncbi.nlm.nih.gov/genome/ - enter species name - click on Genome Assembly and Annotation report - sort table by clicking on the column header Size (Mb) - look at range of sizes in this column.

!!! note "Question"
    Where could you view the output **filename.gfa** and what would it show?

!!! success ""
    ??? "**Answer**"
        This is the assembly graph. You can view it using the tool "Bandage", https://rrwick.github.io/Bandage/, to see how the contigs are connected (including ambiguities).

## Get data

    cd /home/trainee/long_reads/workshop_files
    ls -l

The files we need are:

- **pacbio.fastq.gz</fn>** - the PacBio reads
- **<fn>illumina_R1.fastq.gz</fn>** - the Illumina forward reads
- **<fn>illumina_R2.fastq.gz</fn>** - the Illumina reverse reads

<!--

If you already have the files, skip forward to next section, [Assemble](#assemble).

Otherwise, this section has information about how to find and move the files:

### PacBio files


- Navigate to or create the directory in which you want to work.
- If the files are already on your server, you can symlink by using

```text
ln -s real_file_path [e.g. data/sample_name/pacbio1.fastq.gz] chosen_symlink_name [e.g. pacbio1.fastq.gz]
```

- Alternatively, obtain the input files from elsewhere, e.g. from the BPA portal. (You will need a password.)

- Pacbio files are often stored in the format:
    - **<fn>Sample_name/Cell_name/Analysis_Results/long_file_name_1.fastq.gz</fn>**

- We will use the **<fn>longfilename.subreads.fastq.gz</fn>** files.

- The reads are usually split into three separate files because they are so large.

- Right click on the first **<fn>subreads.fastq.gz</fn>** file and "copy link address".

- In the command line, type:

```text
wget --user username --password password [paste link URL for file]
```
- Repeat for the other two **<fn>subreads.fastq.gz</fn>** files.
- Join the files:
```text
cat pacbio*.fastq.gz > pacbio.fastq.gz
```
- If the files are not gzipped, type:
```text
cat pacbio*.fastq | gzip > pacbio.fastq.gz
```

### Illumina files

- We will also use 2 x Illumina (Miseq) fastq.gz files.
- These are the **<fn>R1.fastq.gz</fn>** and **<fn>R2.fastq.gz</fn>** files.
- Symlink or "wget" these files as described above for PacBio files.
- Shorten the name of each of these files:

```text
mv longfilename_R1.fastq.gz illumina_R1.fastq.gz
mv longfilename_R2.fastq.gz illumina_R2.fastq.gz
```
-->

### NGS Workshop

Some of these tools will take too long to run in this workshop. For these tools, we have pre-computed the output files. In this workshop, we will still enter in the commands and set the tool running, but will sometimes then stop the run and move on to pre-computed output files.

In your directory, along with the PacBio and Illumina files, you may also see folders of pre-computed data.

<!-- list of precomputed data: illuminaR1 and R2, pacbio.fastq.gz, canu_outdir, circlator_outdir, aln.bam, spades_assembly, aln_illumina_pacbio.bam, pilon1.fasta, pilon1.changes -->

### Sample information

The sample used in this tutorial is a gram-positive bacteria called *Staphylococcus aureus* (sample number 25747). This particular sample is from a strain that is resistant to the antibiotic methicillin (a type of penicillin). It is also called MRSA: methicillin-resistant *Staphylococcus aureus*. It was isolated from (human) blood and caused bacteraemia, an infection of the bloodstream.

## Assemble<a name="assemble"></a>

- We will use the assembly software called Canu, <https://github.com/marbl/canu>.
- As we don't have time for Canu to complete, we will look at pre-computed data in the folder **canu_outdir**.

!!! failure "Precomputed section" 
    ```text
    canu -p canu -d canu_outdir_NGS genomeSize=2.8m -pacbio-raw pacbio.fastq.gz    
    ```

    - the first `canu` tells the program to run    
    - `-p canu` names prefix for output files ("canu")    
    - `-d canu_outdir_NGS` names output directory    
    - `genomeSize` only has to be approximate.    
        - e.g. *Staphylococcus aureus*, 2.8m    
        - e.g. *Streptococcus pyogenes*, 1.8m    

    - Canu will correct, trim and assemble the reads.    
    - Various output will be displayed on the screen.    


### Canu output

Move into **<fn>canu_outdir</fn>** and `ls` to see the output files.

```text
cd canu_outdir
ls -l
```

- The **<fn>canu.contigs.fasta</fn>** are the assembled sequences.
- The **<fn>canu.unassembled.fasta</fn>** are the reads that could not be assembled.
- The **<fn>canu.correctedReads.fasta.gz</fn>** are the corrected Pacbio reads that were used in the assembly.
- The **<fn>canu.contigs.gfa</fn>** is the graph of the assembly.
- Display summary information about the contigs: (`infoseq` is a tool from [EMBOSS](http://emboss.sourceforge.net/index.html))

```text
infoseq canu.contigs.fasta
```

!!! note "Question"
    How long is the assembled contig ?

!!! success ""
    ??? "**Answer**"
        tig00000001       2851805
        This looks like a chromosome of approximately 2.8 million bases.

This matches what we would expect for this sample. For other data, Canu may not be able to join all the reads into one contig, so there may be several contigs in the output. Also, the sample may contain some plasmids and these may be found full or partially by Canu as additional contigs.  

!!! hint "Try it later" 
    Change Canu parameters if required: If the assembly is poor with many contigs, re-run Canu with extra sensitivity parameters; e.g.    

    ```text
    canu -p prefix -d outdir corMhapSensitivity=high corMinCoverage=0 genomeSize=2.8m -pacbio-raw pacbio.fastq.gz    
    ```

## Trim and circularise

### Run Circlator

Circlator (<https://github.com/sanger-pathogens/circlator>)
identifies and trims overhangs (on chromosomes and plasmids) and orients the start position at an appropriate gene (e.g. dnaA). It takes in the assembled contigs from Canu, as well as the corrected reads prepared by Canu.

Overhangs are shown in blue:

![circlator](images/circlator_diagram.png)
*Adapted from Figure 1. Hunt et al. Genome Biology 2015*

Run Circlator:

We have already run the command below and now we can look at pre-computed data in the folder **circlator_outdir**    

!!! failure "Pre-computed section" 
    To run the command move back into your main analysis folder.

    ```text
    cd /home/trainee/long_reads/workshop_files
    circlator all --threads 4 --verbose canu_outdir/canu.contigs.fasta canu_outdir/canu.correctedReads.fasta.gz circlator_outdir_NGS    
    ```

    `--threads` is the number of cores <!-- change this to an appropriate number-->    
    `--verbose` prints progress information to the screen    
    `canu_outdir/canu.contigs.fasta` is the file path to the input Canu assembly    
    `canu_outdir/canu.correctedReads.fasta.gz` is the file path to the corrected Pacbio reads - note, fastA not fastQ    
    `circlator_outdir_NGS` is the name of the output directory.    


### Circlator output

Move into the **<fn>circlator_outdir</fn>** directory and `ls` to list files.

```text
cd /home/trainee/long_reads/workshop_files/circlator_outdir
ls -ltr
```


### Questions

!!! note "Question"
    Were all the contigs circularised? Why/why not?

!!! hint ""
    ??? "Hint"
        ```text
        less 04.merge.circularise.log
        ```
        Type "q" to exit.

!!! success ""
    ??? "**Answer**"
        - Yes, the contig was circularised (last column).
        In this example, the contig could be circularized because it contained the entire sequence, with overhangs that were trimmed.

!!! note "Question"
    *Where were the contigs oriented (which gene)?*

!!! hint ""
    ??? "Hint"
        ```text
        less 06.fixstart.log
        ```
        Type "q" to exit.

!!! success ""
    ??? "**Answer**"
        - Look in the "gene_name" column.
        - The contig has been oriented at tr|A0A090N2A8|A0A090N2A8_STAAU, which is another name for dnaA. <!-- (search swissprot - uniprot.org) --> This is typically used as the start of bacterial chromosome sequences.

!!! note "Question"
    *What are the trimmed contig sizes?*

!!! hint ""
    ??? "Hint"
        ```text
        infoseq 06.fixstart.fasta
        ```
!!! success ""
    ??? "**Answer**"
        - tig00000001 2823331 (28564 bases trimmed)
        -This trimmed part is the overlap.

!!! note "Question"
    Circlator can set the start of the sequence at a particular gene. Which gene does it use? Is this appropriate for all contigs?

!!! success ""
    ??? "**Answer**"
        Circlator uses dnaA for the chromosomal contig. For other contigs, it uses a centrally-located gene. However, ideally, plasmids would be oriented on a gene such as a rep gene. It is possible to provide a file to Circlator to do this.

*Re-name the contigs file*:

- The trimmed contigs are in the file called **<fn>06.fixstart.fasta</fn>**.
- Re-name it **<fn>contig1.fasta</fn>**:

```text
cp 06.fixstart.fasta contig1.fasta
```

Open this file in a text editor (e.g. nano: `nano contig1.fasta`) and change the header to ">chromosome".

Move the file back into the main folder (`mv contig1.fasta ../`).

!!! hint "Tips" 
    If all the contigs have not circularised with Circlator, an option is to change the `--b2r_length_cutoff` setting to approximately 2X the average read depth.


## Find smaller plasmids

Pacbio reads are long, and may have been longer than small plasmids. We will look for any small plasmids using the Illumina reads.

This section involves several steps:

1. Use the Canu+Circlator output of a trimmed assembly contig.
2. Map all the Illumina reads against this PacBio-assembled contig.
3. Extract any reads that *didn't* map and assemble them together: this could be a plasmid, or part of a plasmid.
5. Look for overhang: if found, trim.

### Align Illumina reads to the PacBio contig

- Index the contigs file:

```text
bwa index contig1.fasta
```

- Align Illumina reads using using bwa mem:
  We will use the pre-computed file called **aln.bam**.     

!!! failure "Precomputed section" 
    ```text
    bwa mem -t 4 contig1.fasta illumina_R1.fastq.gz illumina_R2.fastq.gz | samtools sort > aln_NGS.bam     
    ```     

    - `bwa mem` is the alignment tool     
    - `-t 4` is the number of cores <!-- choose an appropriate number -->     
    - `contig1.fasta` is the input assembly file     
    - `illumina_R1.fastq.gz illumina_R2.fastq.gz` are the Illumina reads     
    - ` | samtools sort` pipes the output to samtools to sort     
    - `> aln_NGS.bam` sends the alignment to the file **<fn>aln_NGS.bam</fn>**     


### Extract unmapped Illumina reads

- Index the alignment file:

```text
samtools index aln.bam
```

- Extract the fastq files from the bam alignment - those reads that were unmapped to the Pacbio alignment - and save them in various "unmapped" files:

```text
samtools fastq -f 4 -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq aln.bam
```

- `fastq` is a command that coverts a **<fn>.bam</fn>** file into fastq format
- `-f 4` : only output unmapped reads
- `-1` : put R1 reads into a file called **<fn>unmapped.R1.fastq</fn>**
- `-2` : put R2 reads into a file called **<fn>unmapped.R2.fastq</fn>**
- `-s` : put singleton reads into a file called **<fn>unmapped.RS.fastq</fn>**
- `aln.bam` : input alignment file

We now have three files of the unampped reads: **<fn> unmapped.R1.fastq</fn>**, **<fn> unmapped.R2.fastq</fn>**, **<fn> unmapped.RS.fastq</fn>**.

### Assemble the unmapped reads

- Assemble with Spades (<http://cab.spbu.ru/software/spades/>):
- We will use the pre-computed file in the folder **spades_assembly**.

!!! failure "Precomputed section" 
    ```text
    spades.py -1 unmapped.R1.fastq -2 unmapped.R2.fastq -s unmapped.RS.fastq --careful --cov-cutoff auto -o spades_assembly_NGS
    ```

    - `-1` is input file forward
    - `-2` is input file reverse
    - `-s` is unpaired
    - `--careful` minimizes mismatches and short indels
    - `--cov-cutoff auto` computes the coverage threshold (rather than the default setting, "off")
    - `-o` is the output directory


Move into the output directory (**<fn>spades_assembly</fn>**) and look at the contigs:

```text
infoseq contigs.fasta
```

- 78 contigs were assembled, with the max length of 2250 (the first contig).  
- All other nodes are < 650kb so we will disregard as they are unlikely to be plasmids.
- We will extract the first sequence (NODE_1):

```text
samtools faidx contigs.fasta
```

```text
samtools faidx contigs.fasta NODE_1_length_2550_cov_496.613 > contig2.fasta
```

- This is now saved as **<fn>contig2.fasta</fn>**
- Open in nano and change header to ">plasmid".

### Trim the plasmid

To trim any overhang on this plasmid, we will blast the start of contig2 against itself.

- Take the start of the contig:

```text
head -n 10 contig2.fasta > contig2.fa.head
```

- We want to see if it matches the end (overhang).
- Format the assembly file for blast:

```text
makeblastdb -in contig2.fasta -dbtype nucl
```

- Blast the start of the assembly (.head file) against all of the assembly:

```text
blastn -query contig2.fa.head -db contig2.fasta -evalue 1e-3 -dust no -out contig2.bls
```

- Look at **<fn>contig2.bls</fn>** to see hits:

```text
less contig2.bls
```

- The first hit is at start, as expected.
- The second hit is at 2474 all the way to the end - 2550.
- This is the overhang.
- Trim to position 2473.
- Type 'q' to exit.
- Index the plasmid.fa file:

```text
samtools faidx contig2.fasta
```

- Trim

```text
samtools faidx contig2.fasta plasmid:1-2473 > plasmid.fa.trimmed
```

- `plasmid` is the name of the contig, and we want the sequence from 1-2473.

- Open this file in nano (`nano plasmid.fa.trimmed`) and change the header to ">plasmid", save.
- (Use the side scroll bar to see the top of the file.)
- We now have a trimmed plasmid.
- Move file back into main folder:

```text
cp plasmid.fa.trimmed ../
```

- Move into the main folder.

### Plasmid contig orientation

The bacterial chromosome was oriented at the gene dnaA. Plasmids are often oriented at the replication gene, but this is highly variable and there is no established convention. Here we will orient the plasmid at a gene found by Prodigal, in Circlator:

```text
circlator fixstart plasmid.fa.trimmed plasmid_fixstart
```

- `fixstart` is an option in Circlator just to orient a sequence.
- `plasmid.fa.trimmed` is our small plasmid.
- `plasmid_fixstart` is the prefix for the output files.

View the output:

```text
less plasmid_fixstart.log
```

- The plasmid has been oriented at a gene predicted by Prodigal, and the break-point is at position 1200.
- Change the file name:

```text
cp plasmid_fixstart.fasta contig2.fasta
```

<!-- note: annotated with prokka. plasmid only has two proteins. ermC, and a hypothetical protein. protein blast genbank: matches a staph replication and maintenance protein. -->

### Collect contigs

```text
cat contig1.fasta contig2.fasta > genome.fasta
```

- See the contigs and sizes:

```text
infoseq genome.fasta
```

- chromosome: 2823331
- plasmid: 2473

### Questions

!!! note "Question"
    Why is this section so complicated?

!!! success ""
    ??? "**Answer**"
        Finding small plasmids is difficult for many reasons! This paper has a nice summary: On the (im)possibility to reconstruct plasmids from whole genome short-read sequencing data. doi: https://doi.org/10.1101/086744

!!! note "Question"
    Why can PacBio sequencing miss small plasmids?

!!! success ""
    ??? "**Answer**"
        Library prep size selection

!!! note "Question"
    We extract unmapped Illumina reads and assemble these to find small plasmids. What could they be missing?

!!! success ""
    ??? "**Answer**"  
        Repeats that have mapped to the PacBio assembly.

!!! note "Question"
    How do you find a plasmid in a Bandage graph?

!!! success ""
    ??? "**Answer**"
        It is probably circular, matches the size of a known plasmid, and has a rep gene.

!!! note "Question"
    Are there easier ways to find plasmids?

!!! success ""
    ??? "**Answer**"
        Possibly. One option is the program called Unicycler which may automate many of these steps. https://github.com/rrwick/Unicycler

## Correct

We will correct the Pacbio assembly with Illumina reads, using the tool Pilon (<https://github.com/broadinstitute/pilon/wiki>).

### Make an alignment file

- Align the Illumina reads (R1 and R2) to the draft PacBio assembly, e.g. **<fn>genome.fasta</fn>**:
- We will use the pre-computed file called **aln_illumina_pacbio.bam**.


!!! failure "Precomputed section" 
    ```text
    bwa index genome.fasta    
    bwa mem -t 4 genome.fasta illumina_R1.fastq.gz illumina_R2.fastq.gz | samtools sort > aln_illumina_pacbio_NGS.bam    
    ```

    - `-t` is the number of cores <!-- set this to an appropriate number. (To find out how many you have, `grep -c processor /proc/cpuinfo`) -->    



- Index the files:

```text
samtools index aln_illumina_pacbio.bam
samtools faidx genome.fasta
```

- Now we have an alignment file to use in Pilon: **<fn>aln_illumina_pacbio.bam</fn>**

### Run Pilon
Pilon is a software tool which can be used to:

- Automatically improve draft assemblies    
- Find variation among strains, including large event detection    

- We will use the pre-computed files called with the prefixes **pilon1**._

!!! failure "Precomputed section" 
    ```text
    pilon --genome genome.fasta --frags aln_illumina_pacbio.bam --output pilon1_NGS --fix all --mindepth 0.5 --changes --verbose --threads 4
    ```

    - `--genome` is the name of the input assembly to be corrected     
    - `--frags` is the alignment of the reads against the assembly     
    - `--output` is the name of the output prefix     
    - `--fix` is an option for types of corrections     
    - `--mindepth` gives a minimum read depth to use     
    - `--changes` produces an output file of the changes made     
    - `--verbose` prints information to the screen during the run     
    - `--threads`: the number of cores     



Look at the changes file:

```text
less pilon1.changes
```

*Example:*

![pilon](images/pilon.png)


Look at the details of the fasta file:

```text
infoseq pilon1.fasta
```

- chromosome - 2823340 (net +9 bases)
- plasmid - 2473 (no change)


**Option:**

If there are many changes, run Pilon again, using the **<fn>pilon1.fasta</fn>** file as the input assembly, and the Illumina reads to correct.


### Genome output

- Change the file name:

```text
cp pilon1.fasta assembly.fasta
```

- We now have the corrected genome assembly of *Staphylococcus aureus* in .fasta format, containing a chromosome and a small plasmid.  

### Questions

!!! note "Question"
    Why don't we correct earlier in the assembly process?

!!! success ""
    ??? "**Answer**"
        We need to circularise the contigs and trim overhangs first.

!!! note "Question"
    Why can we use some reads (Illumina) to correct other reads (PacBio) ?

!!! success ""
    ??? "**Answer**"
        Illumina reads have higher accuracy

!!! note "Question"
    Could we just use PacBio reads to assemble the genome?

!!! success ""
    ??? "**Answer**"
        Yes, if accuracy adequate.

<!-- ## Short-read assembly: a comparison

So far, we have assembled the long PacBio reads into one contig (the chromosome) and found an additional plasmid in the Illumina short reads.

If we only had Illumina reads, we could also assemble these using the tool Spades.

You can try this here or try it later on your own data.

### Get data

We will use the same Illumina data as we used above:

- **<fn>illumina_R1.fastq.gz</fn>**: the Illumina forward reads
- **<fn>illumina_R2.fastq.gz</fn>**: the Illumina reverse reads

### Assemble

Run Spades:

```text
spades.py -1 illumina_R1.fastq.gz -2 illumina_R2.fastq.gz --careful --cov-cutoff auto -o spades_assembly_all_illumina
```

- `-1` is input file of forward reads
- `-2` is input file of reverse reads
- `--careful` minimizes mismatches and short indels
- `--cov-cutoff auto` computes the coverage threshold (rather than the default setting, "off")
- `-o` is the output directory

### Results

Move into the output directory and look at the contigs:

```text
infoseq contigs.fasta
```

### Questions

How many contigs were found by Spades?

- many

How does this compare to the number of contigs found by assembling the long read data with Canu?

- many more.

Does it matter that an assembly is in many contigs?

- Yes

  - broken genes => missing/incorrect annotations
  - less information about structure: e.g. number of plasmids

- No

  - Many or all genes may still be annotated
  - Gene location is useful (e.g. chromosome, plasmid1) but not always essential (e.g. presence/absence of particular resistance genes)

How can we get more information about the assembly from Spades?

- Look at the assembly graph **<fn>assembly_graph.fastg</fn>**, e.g. in the program Bandage. This shows how contigs are related, albeit with ambiguity in some places.

-->


## Next

### Further analyses

- Annotate genomes, e.g. with Prokka, https://github.com/tseemann/prokka
- Comparative genomics, e.g. with Roary, https://sanger-pathogens.github.io/Roary/

### Links

<!--- [Details of bas.h5 files](https://s3.amazonaws.com/files.pacb.com/software/instrument/2.0.0/bas.h5+Reference+Guide.pdf)
-->

- Canu [manual](http://canu.readthedocs.io/en/stable/quick-start.html) and [gitub repository](https://github.com/marbl/canu)
- Circlator [article](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0849-0) and [github repository](http://sanger-pathogens.github.io/circlator/)
- Pilon [article](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112963) and [github repository](https://github.com/broadinstitute/pilon/wiki)
- Notes on [finishing](https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Finishing-Bacterial-Genomes) and [evaluating](https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/Evaluating-Assemblies) assemblies.
