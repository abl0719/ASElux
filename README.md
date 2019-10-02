Copyright ©2017. University of California, Los Angeles

Created by Zong Miao

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

--------------------------------------------------------------------------------
ASElux: the ultrafast ASE alignment tool

ASElux is a c++ program and it requires a gcc version >= 5. It is designed to 
count the allelic reads overlapping known SNPs from RNA-seq data. To use ASElux,
first, build the static index with a fasta file and a gtf file. The fasta file 
is the sequence of the reference genome in fasta format. The gtf file is the 
gene annotation of the reference genome. In the second step, a vcf file and a 
fasta/fastq file are required to count the alleic reads. The vcf file can be 
obtained from SNP array or genome sequencing. We do not recommend to retrieve 
the vcf file directly from RNA-seq data. Although ASElux works for various 
length of reads, all the reads have to be in the same length in each run.

The output file of ASElux has four colunms: SNP ID, reference allele count,
alternative allele count,and proportion of the reference allele. The SNPs are 
not sorted in any specific order.

Build the program:
A pre-compiled binary file for linux 64bit system is attached and named ASElux.

Compiling from source files requires a gcc version not older than gcc-5. Please 
replace the “gcc” in the Makefile with your gcc directory.


Demo:
There is a random generated demo data set in the ./demo folder.
The demo can be used by running the bash.sh or following the instructions below.

Use:
ASElux [run mode] [--parameters]

build static index:
ASElux build --gtf file --ref file --out filePrefix

alignment:
ASElux align --fa/fq --se/pe --readLen n --index filePrefix --vcf file 
    --seqFiles file1_r1 file1_r2 file2_r1 file2_r2 ... --out fileFrefix

parameters:
fa/fq
    fa      files of reads are fasta format
    fq      files of reads are fastq format

se/pe
    se      use single end reads
    pe      use paired end reads

run mode:
    align   run alignment
    build   build static index

--readLen
    n       the read length of fa/fq files, it has to be indicated by user.
            Although ASElux dose not has a limit set for the read length, reads 
            that are too long (thousands of bp) may cost too much memory and
            time to align.

—-mis
    n       the number of mismatches allowed in each read (default = 2). Since 
            the known SNPs won't be counted as mismatches during the alignment, 
            we do not recommend to set a very high mismatch threshold for a 
            better accuracy.

--index file prefix
    file    the name of the static index

--vcf
    file    the file of vcf format SNPs. Please make sure that vcf file fits the v4.2 standards and contains the GT field.

--seqFiles
    file    the files of reads seperated by space. 
            all files should belongs to the same individual.

--gtf
    file    the annotatino file of the reference genome in gtf format

--ref
    file    the fasta file of the reference genome

--out
    file    the prefix of output file
    
--nthread
    n       the number of threads used during alignment (default = 1).

Update: 2018.04.17
In the default alignment mode, the read counts are SNP based which means that a read
can be counted multiple times if there are more than one SNP existing in the
read. We now offers a new option to count the reads toward only one random SNP in the
read.

New parameter:

--count_once
    for the reads that countain multiple SNPs, ASElux randomly pick one SNP in 
    the read and count towards that SNP to avoid the double counting problem.
