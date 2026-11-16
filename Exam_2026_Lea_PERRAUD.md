---
title: "Exam 2026 - Bioinformatics applied to genomics"
author: Léa PERRAUD
output: html_document
date: "2026-11-13"
---



## 0 - Data
```Bash
cd $HOME

# The files .fastq were downloaded, unzipped and then charged into data/
mkdir data
# We will retrieve them from data folder for the upcoming analysis
# We then create the work directory, where the results will be stored
mkdir work
```

## Exercice 1: Take inventory of the dataset

1) Within data/, let's look at the first record of popA_ind01.fastq with head and identify the four lines using this command: `head -4 popA_ind01.fastq`:
- First line is the header: it starts with @ and contains the read name, the sequencing run information and the read length
- Second line is the nucleotide sequence of the read, which begins with TGCA, the restriction site motif
- Third line is a separator that starts with + and repeats the read name
- Fourth line holds the quality scores, one character per base, so it has the same length as the sequence. Most characters are I (very high quality 40)

2) Let's count the lines of popA_ind01.fastq and compute the number of reads. 
I stored it in a variable n and displayed it with `echo.wc -l popA_ind01.fastq`
```Bash
n=$(($(wc -l < popA_ind01.fastq)/4))
echo "$n"
```
From it, we obtained 300 reads.

3) To answer this question, I used nano, a text editor that runs in the terminal.
This .sh file was created in the work directory (since the data/ directory should have data only)
```bash
cd work/
nano ex1_count_reads.sh 
```
Once done, a text editor opened in the terminal. It is time to type the following loop:
```bash
#!/bin/bash # this line allows to specify the used language, which is bash 
for f in "$@" ; do
    pop=$(echo "$f" | cut -d '_' -f 1)
    n=$(( $(wc -l < "$f") / 4 ))
    echo -e "$f\t$pop\t$n"
done
```
This script takes any number of files as arguments and prints one line per file, separated by tabs.
Once it is done, I close ex1_count_reads.sh and check if the loop has been successfully included: `cat ex1_count_reads.sh # yes !`
Now, let's make it executable `chmod a+x ex1_count_reads.sh`

4) Before anything, let's test on one file first, then I will loop over all of them (as suggested in 0 - Data):
```Bash
./ex1_count_reads.sh ~/data/popA_ind01.fastq # it works indeed ! We obtained file name in column 1, population name in column 2 and the nbr of reads in column 3
```
The screen displays 300 reads !

I ran it on the files of population A only:
```Bash
./ex1_count_reads.sh ~/data/popA_*.fastq
```
We obtain 300, 250 and 200 reads for ind01, ind02 and ind03 respectively, within population A.

Then run it on all files and store them in the work/ directory, as asked in the 0 - Data section:
```Bash
./ex1_count_reads.sh ~/data/pop*.fastq >> ../work/ex1_results.txt # it works ! We typed a star after "pop"" because every filename start with pop, so by doing "pop*" we make sure to include all files
```

5a) How many reads in total for population A? 
```Bash
echo $(($(cat ~/data/popA_*.fastq | wc -l)/4)) # 750 reads
```
We obtain 750 reads for population A.

5b) How many reads in total for population B ? 
```Bash
echo $(($(cat ~/data/popB_*.fastq | wc -l)/4)) # 900
```
We obtain 900 reads for population B.

In total we have:
```Bash
echo $(($(cat pop*.fastq | wc -l)/4)) 
```
We obtain a total of 1650 reads for both populations cumulated.

6) Let's compare the following run with my count:
```Bash
grep -c "^@" ~/data/popA_ind01.fastq # 307 reads
```

This count has more reads (307) than in popA_ind01 (300) because the symbol @ can be read as a quality score character, and some quality scores start with @. That means 7 lines were counted as if they were read names

Let's try it differently:
```Bash
grep -c "^@SRR" ~/data/popA_ind01.fastq # 300 reads
```

We conclude that counting reads with 'grep' is not so precise. Indeed, ^@ also matches quality lines that start with @, and ^@SRR only works because the read names start with SRR.
The safest way is to count the lines and divide by 4, since every FASTQ read is exactly 4 lines.


## Exercice 2: Which individuals look unusable ?

1) For popA_ind01.fastq, I extracted the sequence lines and counted how many contain no N:
```Bash
awk 'NR%4==2' ~/data/popA_ind01.fastq | grep -vc "N" # 278 reads
```

2) Let's compute the percentage of reads without any N:
```Bash
echo $(( 100 * $(awk 'NR%4==2' ~/data/popA_ind01.fastq | grep -vc "N") / ($(wc -l < ~/data/popA_ind01.fastq) / 4) )) # 100* (NpopA/Ntot)
```
We have 92% reads without any N !

3) The ex2_quality_N.sh file was created in the work/ directory as followed:
```Bash
cd work/
nano ex2_quality_N.sh
```

We then copy the following code in the ex2_cutsite.sh
```Bash
#!/bin/bash # this line allows to specify the used language, which is bash
printf "file\tn_reads\tn_reads_without_N\tpercent\n"
for f in "$@" ; do
    n=$(( $(wc -l < "$f") / 4 ))
    sans_N=$(awk 'NR%4==2' "$f" | grep -vc "N")
    printf "%s\t%s\t%s\t%s\n" "$f" "$n" "$sans_N" $(( 100 * sans_N / n ))
done
```

I make sure it is executable !
```Bash
chmod a+x ex2_quality_N.sh   
```

4) Before anything, let's test on one file first (as suggested in 0 - Data):
`./ex2_quality_N.sh ~/data/popA_ind01.fastq`
Let's now see the results for every population !
`./ex2_quality_N.sh ~/data/pop*.fastq` 

The percentage is greater or 90% except in popB_ind03.fastq.
Indeed, popB_ind3 that has 25% reads, which means it needs to be removed from a population genetic analysis.

Store the results in work/ (as suggested in the 0 - Data section):
`./ex2_quality_N.sh ~/data/pop*.fastq >> ../work/ex2_results.txt`


## Exercice 3: Check the restriction site

Context: in RAD-seq, DNA is cut by a restriction enzyme, so nearly all reads should start 
with the same short motif. In the Pedicularis data of the course this motif is TGCA. 

Reads that do not start with it come from something else (contamination, adapters, a failed digest).
1) I would like to count the reads of popA_ind01.fastq whose sequence starts with TGCA
```Bash
grep -c "^TGCA" ~/data/popA_ind01.fastq
```
We obtain 284 reads in individual 01 of population A.

2) To answer this question, I created an ex3_cutsite.sh using the nano command as before (I think it's a nice way of doing things).
As always, I stored the filed in the work/ directory:
```Bash
cd work/
nano ex3_cutsite.sh
```

And then copied the following code in the ex3_cutsite.sh
```Bash
#!/bin/bash # this line allows to specify the used language, which is bash
motif="$1"
shift
printf "file\tmotif\tn_reads_with_motif\tpercent"
for f in "$@" ; do
    n=$(( $(wc -l < "$f") / 4 ))
    avec=$(grep -c "^$motif" "$f")
    printf "%s\t%s\t%d\t%d%%\n" "$f" "$motif" "$avec" $(( 100 * avec / n ))
done 
```
Here, ex3_cutsite.sh takes the motif as first argument and the files as the following arguments.

Of course, we make it executable:
```Bash
chmod a+x ex3_cutsite.sh
```

3) The output contains 4 columns: first is for the file, second is for the motif, third is for n_reads_with_motifs, fourth is the pourcentage

4) Before anything, let's test on one file first, then loop over all of them (as suggested in 0 - Data):
```Bash
./ex3_cutsite.sh TGCA ~/data/popA_ind01.fastq # it works indeed ! We obtained file name in column 1, population name in column 2 and the nbr of reads in column 3
```
Let's now check all populations at once:
```Bash
./ex3_cutsite.sh TGCA ~/data/pop*.fastq
```
We then store the results in work/
```Bash
./ex3_cutsite.sh TGCA ~/data/pop*.fastq >> ../work/ex3_results.txt
```
So, Individual 03 from population B behaves differently as only 23 % of their reads starts with TGCA, against 93 - 97 % for the five others.
It also only has 52 reads with motif, against at least 194 for the other populations.
Two possible explanations:
- Incomplete digestion, indeed, if the enzyme cut in the wrong way, the fragments did not start at the restriction site and therefore, many reads may come from random cuts
- Contamination or bad DNA quality: some reads may come from another organism (as mentioned in class with Amelia) or degraded DNA, that does not contain TGCA sites at the start of the read


## Exercice 4: A small test dataset, and archiving

Context: you never test a new pipeline on the full dataset. You build a small subsample first, and you archive and compress what you send to collaborators

1) All right, I guess it's time to create the folder subsample. Since it was not specified in
the instructions where to create the folder, I decided it would be inside the work/ directory:
```Bash
cd work/
mkdir subsample
```

2) I now want to write ex4_subsample.sh taking the number of reads to keep as $1 and the files as the other arguments. 
```Bash
cd work/
nano ex4_subsample.sh
```

Inside of it, I wrote the following code:
```Bash
#!/bin/bash # this line allows to specify the used language, which is bash
n=$1
shift
lines=$((4*n)) # a read is 4 lines

for f in "$@" ; do
    name=$(basename "$f")
    head -n "$lines" "$f" > subsample/"$name"
done
```

We make it executable (as suggested in the 0 - Data section):
`chmod a+x ex4_subsample.sh`

Check all populations at once (the populations come from data/):
`./ex4_subsample.sh 100 ~/data/pop*.fastq`
I tested with the 100 first reads.

3) Let's check the results in the work/ directory using `wc -l subsample/*`: 

- 400 subsample/popA_ind01.fastq
- 400 subsample/popA_ind02.fastq
- 400 subsample/popA_ind03.fastq
- 400 subsample/popB_ind01.fastq
- 400 subsample/popB_ind02.fastq
- 400 subsample/popB_ind03.fastq

With 100 reads, it matches 400 lines/file. Good !

Let's store the results in the work/ directory:
`wc -l subsample/* > ~/work/ex4_results.txt`

4) It is now time to make a compressed archive subsample.tar.gz of the folder with tar, as seen in class
```Bash
cd work/
tar -czvf subsample.tar.gz subsample/
```

5) In work/, I compared `du -h subsample` with `ls -ls subsample.tar.gz`:
- subsample is 148 K
- subsample.tar.gz is 20 K
Therefore, the archive is about 7.4 times smaller because gzip compresses FASTQ text well.

6) Lastly, I compressed a copy of the `popA_ind01.fastq` file with gzip and compared the sizes with `ls -lh`. 
```Bash
cp ~/data/popA_ind01.fastq ~/work/
```
Now that I am in the work/ directory, I type:
```Bash
ls -lh popA_ind01.fastq
```
The resulting size is 62 K

I am still in the work/ directory,
```Bash
gzip popA_ind01.fastq
ls -lh popA_ind01.fastq.gz
```
The resulting size is 11 K

Therefore, the compression ratio is 5.6. For a project with 500 individuals, fastq files are a few GB,so it
is important to compress them, otherwise we could loose a large amount of disk space. Besides, compressing
files is a good way of making transfers and backups faster, and it also is less expensive.