Key Learning Outcomes
---------------------

After completing this practical the trainee should be able to:

-   Familiarise yourself with the command line environment on a Linux
    operating system.

-   Run some basic linux system and file operation commands

-   Navigation of biological data files structure and manipulation

## Resources

### Tools

* Basic Linux system commands on an Ubuntu OS.

* Basic file operation commands

### Links

* [Software Carpentry](http://software-carpentry.org)

* [Example 1000Genome Project data](http://www.1000genomes.org)

## Author Information

*Primary Author(s):*
    Matt Field matt.field@anu.edu.au     
 

Shell Exercise
--------------

Let’s try out your new shell skills on some real data.

The file `1000gp.vcf` is a small sample (1%) of a very large text file
containing human genetics data. Specifically, it describes genetic
variation in three African individuals sequenced as part of the 1000
Genomes Project (http://www.1000genomes.org). The ’vcf’ extension lets
us know that it’s in a specific text format, namely ’Variant Call
Format’. The file starts with a bunch of comment lines (they start with
’\#’ or ’\#\#’), and then a large number of data lines. This VCF file
lists the differences between the three African individuals and a
standard ’individual’ called the reference (actually based upon a few
different people). Each line in the file corresponds to a difference.
The line tells us the position of the difference (chromosome and
position), the genetic sequence in the reference, and the corresponding
sequence in each of the three Africans. Before we start processing the
file, let’s get a high-level view of the file that we’re about to work
with.

Open the Terminal and go to the directory where the data are stored:
```
cd /home/trainee/cli
ls
pwd
ls -lh 1000gp.vcf
wc -l 1000gp.vcf
```

!!! note "Question"
    What is the file size (in kilo-bytes), and how many lines are in the file?.
!!! hint ""
    ??? "Hint"
        Hint: `man ls`, `man wc`

!!! success ""
    ??? "Answer"
        3.6M

        45034 lines

Because this file is so large, you’re going to almost always want to
pipe ('|') the result of any command to less (a simple text viewer, type
‘`q`’ to exit) or head (to print the first 10 lines) so that you don’t
accidentally print 45,000 lines to the screen.

Let’s start by printing the first 5 lines to see what it looks like.
```
head -5 1000gp.vcf
```

That isn’t very interesting; it’s just a bunch of the comments at the
beginning of the file (they all start with ’\#’)!

Print the first 20 lines to see more of the file.
```
head -20 1000gp.vcf
```

Okay, so now we can see the basic structure of the file. A few comment
lines that start with ’\#’ or ’\#\#’ and then a bunch of lines of data
that contain all the data and are pretty hard to understand. Each line
of data contains the same number of fields, and all fields are separated
with TABs. These fields are:

1.  the chromosome (which volume the difference is in)

2.  the position (which character in the volume the difference starts
    at)

3.  the ID of the difference

4.  the sequence in the reference human(s)

The rest of the columns tell us, in a rather complex way, a bunch of
additional information about that position, including: the predicted
sequence for each of the three Africans and how confident the scientists
are that these sequences are correct.

To start analyzing the actual data, we have to remove the header.

!!! note "Question"
    How can we print the first 10 non-header lines (those that don’t start
    with a ’\#’)?

!!! hint ""
    ??? "Hint"
        Hint: `man grep` (remember to use pipes '|')

!!! success ""
    ??? "Answer"
        grep -v "^\#" 1000gp.vcf | head  

!!! note "Question"
    How many lines of data are in the file (rather than counting the number
    of header lines and subtracting, try just counting the number of data
    lines)?

!!! success ""
    ??? "Answer"
        grep -v "^\#" 1000gp.vcf | wc -l (should print 45024)

Where these differences are located can be important. If all the
differences between two encyclopedias were in just the first volume,
that would be interesting. The first field of each data line is the name
of the chromosome that the difference occurs on (which volume we’re on).

!!! note "Question"
    Print the first 10 chromosomes, one per line.
!!! hint ""
    ??? "Hint"
        Hint: `man cut` (remember to remove header lines first)

!!! success ""
    ??? "Answer"
        grep -v "\^\#" 1000gp.vcf | cut -f 1 | head

As you should have observed, the first 10 lines are on numbered
chromosomes. Every normal cell in your body has 23 pairs of chromosomes,
22 pairs of ‘autosomal’ chromosomes (these are numbered 1-22) and a pair
of sex chromosomes (two Xs if you’re female, an X and a Y if you’re
male).

Let’s look at which chromosomes these variations are on.

!!! note "Question"
    Print a list of the chromosomes that are in the file (each chromosome
    name should only be printed once, so you should only print 23 lines).
!!! hint ""
    ??? "Hint"
        Hint: remove all duplicates from your previous answer (`man sort`)

!!! success ""
    ??? "Answer"
        grep -v "\^\#" 1000gp.vcf | cut -f 1 | sort -u

Rather than using `sort` to print unique results, a common pipeline is
to first sort and then pipe to another UNIX command, `uniq`. The `uniq`
command takes sorted input and prints only unique lines, but it provides
more flexibility than just using sort by itself. Keep in mind, if the
input isn’t sorted, `uniq` won’t work properly.

!!! note "Question"
    Using `sort` and `uniq`, print the number of times each chromosome
    occurs in the file.
!!! hint ""
    ??? "Hint"
        Hint: `man uniq`

!!! success ""
    ??? "Answer"
        grep -v "\^\#" 1000gp.vcf | cut -f 1 | sort | uniq -c

!!! note "Question"
    Add to your previous solution to list the chromosomes from most
    frequently observed to least frequently observed.
!!! hint ""
    ??? "Hint"
        Hint: Make sure you’re sorting in descending order. By default, sort
        sorts in ascending order.

!!! success ""
    ??? "Answer"
        grep -v "\^\#" 1000gp.vcf | cut -f 1 | sort | uniq -c | sort -n -r

This is great, but biologists might also like to see the chromosomes
ordered by their number (not dictionary order), since different
chromosomes have different attributes and this ordering allows them to
find a specific chromosome more easily.

!!! note "Question"
    Sort the previous output by chromosome number
!!! hint ""
    ??? "Hint"
        Hint: A lot of the power of sort comes from the fact that you can
        specify which fields to sort on, and the order in which to sort them. In
        this case you only need to sort on one field.

!!! success ""
    ??? "Answer"
        grep -v "\^\#" 1000gp.vcf | cut -f 1 | sort | uniq -c | sort -k 2n
