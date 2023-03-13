#! /usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd;
use FindBin qw($Bin);
use Data::Dumper;
my ($file,$Bin_dir,$gene_dir);
my ($outdir,$identity,$length,$mismatch,$gap,$prefix) = ("./",100,500,0,0,"All");
GetOptions (
	'file:s'	=>	\$file,
	'bin:s'	=>	\$Bin_dir,
	'genedir:s'	=>	\$gene_dir,
	'prefix:s'	=>	\$prefix,
	'outdir:s'	=>	\$outdir,
	'id:f'	=>	\$identity,
	'len:i'	=>	\$length,	
	'mis:i'	=>	\$mismatch,
	'gap:i'	=>	\$gap
);
sub usage {
	die "perl $0
 -file	[s]	MetaCHIP *_detected_HGTs.txt
 -bin	[s]	Bins directory
 -genedir[s]	Gene prediction dir
 -prefix[s]	Output prefix [default All]
 -outdir[s]	Output directory [default ./]
 -id	[f]	Blat identity [default 100]
 -len	[i]	Blat aligment length [default 500]
 -mis	[i]	Blat aligment fregment miss number [default 0]
 -gap	[i]	Blat aligment fregment gap number [default 0]\n";
}
&usage and exit unless (defined $file && defined $Bin_dir);
my $pwd = getcwd;
$outdir = ($outdir eq "." || $outdir eq "./") ? $pwd : ($outdir =~ /^\//) ? $outdir : $pwd."/".$outdir;
system ("mkdir -m 755 -p $outdir") unless (-d $outdir);
my $blast_dir = "$outdir/$prefix.Blat_result_add";
system ("mkdir -m 755 -p $index_dir") unless (-d $index_dir);
system ("mkdir -m 755 -p $blast_dir") unless (-d $blast_dir);
open IN, "$file";
open OUT, ">$outdir/$prefix.MetaCHIP_Blat_Merge_BM.xls";
my $shell = "";
while (my $line = <IN>)
{
	chomp ($line);
	next if ($line =~ /^Gene/);
	next if ($line =~ /^#/);
	my @hgt = split /\t/, $line;
	my @query = split /\_/,$hgt[0];
	my @subject = split /_/,$hgt[1];
	my $query_name = join "_", @query[0..($#query-1)];
	my $subject_name = join "_", @subject[0..($#subject-1)];
	my $query_gene_info = &get_gene_info("$gene_dir/$query_name.gbk");
	my $subject_gene_info = &get_gene_info("$gene_dir/$subject_name.gbk");
	my %query_gene_info = %$query_gene_info;
	my %subject_gene_info = %$subject_gene_info;
	#=== run genome vs genome blast ===
	system ("ln -sf $Bin_dir/$subject_name.fa $index_dir");
	system ("$Bin/blat $index_dir/$subject_name.fa $Bin_dir/$query_name.fa -out=blast8 $blast_dir/$query_name\_$subject_name\_blat_result.xls");
	#=== get genome coverage ===
	`awk '\$4>200 && \$5<50' $blast_dir/$query_name\_$subject_name\_blat_result.xls > $blast_dir/$query_name\_$subject_name\_blat_result_filter.xls`;
	system ("$Bin/soap.coverage -cvg -i $blast_dir/$query_name\_$subject_name\_blat_result_filter.xls -refsingle $Bin_dir/$query_name.fa -o $blast_dir/$query_name\_$subject_name\_query_coverage.xls -m8query");
	system ("$Bin/soap.coverage -cvg -i $blast_dir/$query_name\_$subject_name\_blat_result_filter.xls -refsingle $Bin_dir/$subject_name.fa -o $blast_dir/$query_name\_$subject_name\_subject_coverage.xls -m8subject");
	chomp(my $query_coverage = (split /\:/, `grep -E "^Percentage" $blast_dir/$query_name\_$subject_name\_query_coverage.xls`)[1]);
	chomp(my $subject_coverage = (split /\:/, `grep -E "^Percentage" $blast_dir/$query_name\_$subject_name\_subject_coverage.xls`)[1]);	
	#=== get HGT gene position ===
	my ($gene_q_start,$gene_q_end) = split /\.\./,$query_gene_info{$hgt[0]}{'position'};
	my ($gene_s_start,$gene_s_end) = split /\.\./,$subject_gene_info{$hgt[1]}{'position'};
	my $gene_length = abs ($gene_q_start - $gene_q_end + 1);
	#=== get blast result which match to HGT result ===
	open BLAST, "$blast_dir/$query_name\_$subject_name\_blat_result.xls";
	my $Blast_result = "";
	while (<BLAST>)
	{
		chomp;
		my @array = split /\t/;
		next if ($array[2] < $identity || $array[3] < $length); 
		next if ($array[4] > $mismatch || $array[5] > $gap);
		my ($q_start,$q_end,$s_start,$s_end);
		if($array[6] < $array[7])
		{
			$q_start = $array[6];
			$q_end = $array[7];
		} else {
			$q_start = $array[7];
			$q_end = $array[6];
		}
		if ($array[8] < $array[9])
		{
			$s_start = $array[8];
			$s_end = $array[9];
		} else {
			$s_start = $array[9];
			$s_end = $array[8]
		}
		if ($query_gene_info{$hgt[0]}{'contig'} eq $array[0] && $subject_gene_info{$hgt[1]}{'contig'} eq $array[1])
		{
			if ($q_start <= $gene_q_start && $q_end >= $gene_q_end && $s_start <= $gene_s_start && $s_end >= $gene_s_end) 
			{	
				$Blast_result = "$gene_length\t$gene_q_start..$gene_q_end\t$gene_s_start..$gene_s_end\t$array[0]\t$array[1]\t$query_coverage\t$subject_coverage\t$array[2]\t$array[3]\t$q_start..$q_end\t$s_start..$s_end";
				print OUT "$line\t$Blast_result\n";
			}
		}
	}
	if ($Blast_result eq "")
	{
		$Blast_result = "$gene_length\t$gene_q_start..$gene_q_end\t$gene_s_start..$gene_s_end\t$query_gene_info{$hgt[0]}{'contig'}\t$subject_gene_info{$hgt[1]}{'contig'}\t$query_coverage\t$subject_coverage";
		print OUT "$line\t$Blast_result\n";
	}
}
#================ subroutine ===================
sub get_gene_info {
	my ($file) = @_;
	open GENE, "$file";
	$/="LOCUS";<GENE>;
	my %Gene;
	while (<GENE>)
	{
		chomp;
		my @line = split /\n/;
		my $contig = (split/\s+/,$line[0])[1];
		for (my $i=1; $i<=$#line; $i++)
		{
			if($line[$i]=~/^\s+CDS\s+/)
			{
				my $loc=$1 if ($line[$i]=~/(\d+\.\.\d+)/);
				my $gene=$1 if ($line[$i+1]=~/locus_tag=\"(\S+)\"/);
				$Gene{$gene}{'contig'} = $contig;
				$Gene{$gene}{'position'} = $loc;
			}
		}
	}
	close GENE;
	$/="\n";
	return (\%Gene);
}
#============ THE END ============
