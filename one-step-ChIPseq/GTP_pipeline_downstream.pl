use warnings;
use strict;

####################################################
#      Written by Jiyang Li (meng_ljy@126.com)     #
####################################################

die "perl $0 [input_text file samplename, sample.sorted.bam, input.bam] [out_dir]\n" unless @ARGV == 2;

my ($in_f, $out_f) = @ARGV;
die "Overlap In-Output...\n" if $in_f eq $out_f;



##samlename	fastq1.gz	fastq2.gz
my (@arr);
##
if ($in_f =~ /\.gz$/){
	open IN, "gzip -dc $in_f |" or die $!;
}else{
	open IN, $in_f or die $!;
}
##main
#
	while(<IN>){
		chomp;
		my @script;
		my @a =split /\t/,$_;
		my $prefix = $a[0];
		my $samplebam = $a[1];
		my $inputbam = $a[2];
		open OT, ">$out_f/$prefix.downstream.sh" or die $!;
		push @arr, "$out_f/$prefix.downstream.sh >$prefix.bamDownstream.log 2>&1";
		&fileexist('bedfiles');
		&fileexist('gokegg');
		&fileexist('heatmapprofile');
		&fileexist('motif');
		&fileexist('superenhancer');
		push @script, &callpeakspe($prefix, $samplebam, $inputbam);
		push @script, &getpos($prefix, 'broadPeak');
		push @script, &getpos($prefix, 'narrowPeak');
		push @script, &chipanno_gokegg_enrigh($prefix, 'broadPeak');
		push @script, &chipanno_gokegg_enrigh($prefix, 'narrowPeak');
		push @script, &computematrix($prefix, 'broadPeak');
		push @script, &computematrix($prefix, 'narrowPeak');
		push @script, &plotheadmap($prefix, 'broadPeak');
		push @script, &plotheadmap($prefix, 'narrowPeak');
		push @script, &plotprofile($prefix, 'broadPeak');
		push @script, &plotprofile($prefix, 'narrowPeak');
		push @script, &findmotif($prefix,'broadPeak');
		push @script, &findmotif($prefix,'narrowPeak');
		push @script, &findsuperenhancer($prefix,'broadPeak', $samplebam, $inputbam);
    	push @script, &findsuperenhancer($prefix,'narrowPeak', $samplebam, $inputbam);
		my $p = join("\n",@script);
		print OT "$p\n";
		close OT;
	}
	close IN;
open OT, ">$out_f/work.sh" or die $!;
my $q = join("\n",@arr);
print OT "$q\n";
close OT;


sub fileexist{
	my ($f) = @_;
	mkdir("$f")unless(-d "$f");
}

##call peak
sub callpeakspe{
	my ($prefix,$samplebam,$inputbam) = @_;
	my $p1 = "/picb/molsysbio/usr/lijiyang/software/install/anaconda3/bin/macs2 callpeak -t $samplebam -c $inputbam -n ./bedfiles/$prefix.compared -g mm --broad --format BAMPE\n";
	my $p2 = "/picb/molsysbio/usr/lijiyang/software/install/anaconda3/bin/macs2 callpeak -t $samplebam -c $inputbam -n ./bedfiles/$prefix.compared -g mm --format BAMPE -B\n";
	my $p3 = "/picb/molsysbio/usr/lijiyang/software/install/anaconda3/bin/bedtools intersect -a ./bedfiles/$prefix.compared_peaks.broadPeak -b /picb/molsysbio/usr/lijiyang/project/00.ShareInfo/01.pipeline/01.ChIP-seq/blacklists/mm10/mm10.blacklist.bed -v > ./bedfiles/$prefix.compared.broadPeak.filtered.bed\n";
	my $p4 = "/picb/molsysbio/usr/lijiyang/software/install/anaconda3/bin/bedtools intersect -a ./bedfiles/$prefix.compared_peaks.narrowPeak -b /picb/molsysbio/usr/lijiyang/project/00.ShareInfo/01.pipeline/01.ChIP-seq/blacklists/mm10/mm10.blacklist.bed -v > ./bedfiles/$prefix.compared.narrowPeak.filtered.bed";
	my $merge = $p1.$p2.$p3.$p4;
	return $merge;
}

##


sub getpos{
	my ($prefix,$type) = @_;
	my $i = "/picb/molsysbio/usr/lijiyang/software/install/anaconda3/bin/bedops -i ./bedfiles/$prefix.compared.$type.filtered.bed ./bedfiles/$prefix.compared.$type.filtered.bed > ./bedfiles/$prefix.compared.$type.filtered.bed.PosOnly";
	return $i;
}


sub chipanno_gokegg_enrigh{
	my ($prefix,$type) = @_;
	my $r = "Rscript /picb/molsysbio/usr/lijiyang/project/00.ShareInfo/01.pipeline/01.ChIP-seq/script/annoChIPseqByChIPseeker.r ./bedfiles/$prefix.compared.$type.filtered.bed.PosOnly  ./gokegg/$prefix";
	return $r;
}

## plot heatmap and profile
sub computematrix{
	 my ($prefix,$type) = @_;
	my $r = "/picb/molsysbio/usr/lijiyang/software/install/anaconda3/bin/computeMatrix reference-point --referencePoint TSS -b 3000 -a 10000 -R ./bedfiles/$prefix.compared.$type.filtered.bed -S bigwig/$prefix.bw -p 10 --skipZeros -o ./bedfiles/matrix_$prefix.$type.gz --outFileSortedRegions ./bedfiles/SortedRegions.$prefix.$type.bed ";
}

sub plotheadmap{
	my ($prefix,$type) = @_;
	my $r = "/picb/molsysbio/usr/lijiyang/software/install/anaconda3/bin/plotHeatmap -m ./bedfiles/matrix_$prefix.$type.gz -out ./heatmapprofile/$prefix.$type.heatmap.pdf --colorMap RdBu";
	return $r; 
}
sub plotprofile{
	my ($prefix,$type) = @_;
	my $r = "/picb/molsysbio/usr/lijiyang/software/install/anaconda3/bin/plotProfile -m ./bedfiles/matrix_$prefix.$type.gz -out ./heatmapprofile/$prefix.$type.profile.pdf --colorMap RdBu";
	return $r;
}

## type= broadPeak or narrowPeak
sub findmotif{
	 my ($prefix,$type) = @_;
	 my $r =  "perl /picb/molsysbio/usr/lijiyang/software/install/anaconda3/bin/findMotifsGenome.pl ./bedfiles/$prefix.compared.$type.filtered.bed mm10 -p 10 motif/$prefix.$type.motif";
	 return $r;
}

##
sub findsuperenhancer{
	my ($prefix,$type,$samplebam,$inputbam) = @_;
#	my $pwd=`pwd`;
#	&fileexist("$pwd/superenhancer/$prefix.$type");
#	my $r1 = 'cd /home/lijiyang/bin/script/young_computation-rose-1a9bb86b5464'."\n";
#	my $r2 = "python ROSE_main.py -g mm10 -i $pwd/bedfiles/$prefix.compared.$type.filtered.bed -r $samplebam -c $inputbam -o $pwd/superenhancer/$prefix.$type";
#	return $r1.$r2;
	my $r = "perl /home/lijiyang/bin/script/young_computation-rose-1a9bb86b5464/ROSEsuperenhancer.pl $prefix $type $samplebam $inputbam";
	return $r;
}	


