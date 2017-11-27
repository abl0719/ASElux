ASElux="../ASElux_1.0"

#build the static index for ASElux
$ASElux build --gtf annotation.gtf --ref genome.fa --out demo

#align the fasta files
$ASElux align --fa --pe --readLen 50 --index demo --vcf snps.vcf --seqFiles read1.fa read2.fa --out demo_rc.txt
