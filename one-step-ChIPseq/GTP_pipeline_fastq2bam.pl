use warnings;
use strict;

####################################################
#      Written by Jiyang Li (meng_ljy@126.com)     #
#      date:									   #
####################################################

die "perl $0 [input_text file info/samplename/fastq1/fastq2] [out_dir]\n" unless @ARGV == 2;

my ($in_f, $out_f) = @ARGV;
die "Overlap In-Output...\n" if $in_f eq $out_f;

##samlename	fastq1.gz	fastq2.gz
my ($prefix, $fq1, $fq2);
my (@arr);
##
if ($in_f =~ /\.gz$/){
	open IN, "gzip -dc $in_f |" or die $!;
}else{
	open IN, $in_f or die $!;
}

	while(<IN>){
		chomp;
		my @script;
		my @a = split /\t/,$_;
		$prefix = $a[0];
		$fq1 = $a[1];
		$fq2 = $a[2];
		open OT, ">$out_f/$prefix.sh" or die $!;
		push @arr, "$out_f/$prefix.sh >$prefix.bam2fastq.log 2>&1";
		&fileexist('beforeTrim_fastqc');
		&fileexist('afterTrim_fastqc');
		&fileexist('sortedBam');
		&fileexist('bigwig');
		push @script, &fastqc($fq1, $prefix, 'beforeTrim_fastqc');
		push @script, &fastqc($fq2, $prefix, 'beforeTrim_fastqc');
		push @script, &trimmomatic($fq1, $fq2, $prefix);
		push @script, &fastqc("output_$prefix.forward_paired.fq.gz", $prefix, 'afterTrim_fastqc');
		push @script, &fastqc("output_$prefix.reverse_paired.fq.gz", $prefix, 'afterTrim_fastqc');
		push @script, &bowtie2("output_$prefix.forward_paired.fq.gz", "output_$prefix.reverse_paired.fq.gz", $prefix);
		push @script, &sam2bam($prefix, 'sortedBam');
		push @script, &bam2bigwig($prefix);
		my $p = join ("\n", @script);
		print OT "$p\n";
		close OT;
	}
	close IN;

open OT, ">$out_f/work.sh" or die $!;
my $q = join("\n",@arr);
print OT "$q\n";
close OT;

##
sub fileexist{
	my ($f) = @_;
	mkdir("$f")unless(-d "$f");
}

##
sub fastqc{
	my ($fq, $prefix,$dir) = @_;
	my $r = "/picb/molsysbio/usr/lijiyang/software/install/anaconda3/bin/fastqc -o $dir $fq";
	return $r;
}

##
sub trimmomatic{
	my ($fq1, $fq2, $prefix) = @_;
	my $r = "/picb/molsysbio/usr/lijiyang/software/install/anaconda3/bin/trimmomatic PE -phred33 $fq1 $fq2 output_$prefix.forward_paired.fq.gz output_$prefix.forward_unpaired.fq.gz output_$prefix.reverse_paired.fq.gz output_$prefix.reverse_unpaired.fq.gz ILLUMINACLIP:/picb/molsysbio/usr/lijiyang/software/install/anaconda3/share/trimmomatic-0.36-5/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";
	return $r;
}

##
sub bowtie2{
	my ($fq1, $fq2, $prefix) = @_;	
	my $r = "/picb/molsysbio/usr/lijiyang/software/install/anaconda3/bin/bowtie2 -p 4 -x /picb/molsysbio/usr/lijiyang/project/00.ShareInfo/refs/mm10/bowtie2_index/mm10 -1 $fq1 -2 $fq2 -S $prefix.sam";
	return $r;
}

##
sub sam2bam{
	my ($prefix,$dir) = @_;
	my $r1 = "/picb/molsysbio/usr/lijiyang/software/install/anaconda3/bin/samtools view -bS $prefix.sam | samtools sort -@ 4 - -T $prefix -o sortedBam/$prefix.sorted.bam\n";
	my $r2 = "/picb/molsysbio/usr/lijiyang/software/install/anaconda3/bin/samtools index sortedBam/$prefix.sorted.bam\n";
	my $r3 = "rm $prefix.sam";
	my $merge = $r1.$r2.$r3;
	return $merge;
}

##
sub bam2bigwig{
	my ($prefix) = @_;
	my $r = "/picb/molsysbio/usr/lijiyang/software/install/anaconda3/bin/bamCoverage -b sortedBam/$prefix.sorted.bam --normalizeUsing RPKM --binSize 30 --smoothLength 300 -p 4 --extendReads 200 -o bigwig/$prefix.bw";
	return $r;
}


