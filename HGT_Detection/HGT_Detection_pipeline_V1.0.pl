#! /usr/bin/perl
=head1 Description: Detect horizontal gene transfer (HGT) using MetaCHIP and Blat.
 More information: https://github.com/songweizhi/MetaCHIP
 Reference: https://doi.org/10.1186/s40168-019-0649-y; https://doi.org/10.1016/j.cell.2021.02.052
=head1 Usage:
 perl HGT_Detection_pipeline_V1.0.pl
 Required:
 -conf	[s]	config file, please refer to corresponding SOP to get instruction
 -list	[s]	reads pathway and bin pathways, format: SampleName bin_folder_pathwas bin_anno_file checkM_file bin_stat_file
 -outdir[s]	output directory [default: .]
 -softpath[s]	software pathway [default: ./bin]

=head1 Example:
 perl HGT_Detection_pipeline_V1.0.pl -conf HGT_Detection.pl.conf -list Assembly.list
=cut

use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use File::Basename;
use lib "$Bin/lib";
use AnaMethod;
use Data::Dumper;
my ($conf,$list);
my ($run,$outdir) = ("Y",".");
GetOptions (
	'conf:s'	=>	\$conf,
	'list:s'	=>	\$list,
	'outdir:s'	=>	\$outdir
);
die `pod2text $0` unless (defined $conf && defined $list);
chomp (my $pwd = `pwd`);
$outdir = ($outdir eq "." || $outdir eq "./") ? $pwd : ($outdir =~ /^\//) ? $outdir : $pwd."/".$outdir;
my ($dir_list, $dir_process, $dir_sh, $dir_result) = ("$outdir/list", "$outdir/process/HGT_Detection", "$outdir/shell/HGT_Detection", "$outdir/upload/HGT_Detection");
for ($outdir, $dir_list, $dir_process, $dir_sh, $dir_result) 
{
	`mkdir -m 755 -p $_` unless (-d $_);
}
$conf = $pwd."/".$conf if ($conf !~ /^\//);
$list = $pwd."/".$list if ($list !~ /^\//);
#-- read config --
my $config = &readconf($conf);
my %config = %$config;
#== read assembly result ==
open READ, "$list";
my %sample;
my $sample_num;
while (<READ>)
{
	chomp;
	next if (/^#/);
	$sample_num ++;
	my @array = split /\s+/;
	$sample{$array[0]}{'bin'} = $array[1];
	$sample{$array[0]}{'anno'} = $array[2];
	$sample{$array[0]}{'checkm'} = $array[3];
	$sample{$array[0]}{'stat'} = $array[4];
}
close READ;
#===================== product shells ===================
open DEP, ">$dir_list/HGT_Detection_dependence.txt";
open LIST, ">$dir_list/HGT_result.list";
open ASS, ">$dir_process/Species_num_for_HGT.xls";
my $shell = "";
foreach my $sample (sort keys %sample)
{
	system ("mkdir -m 755 -p $dir_process/$sample") unless (-d "$dir_process/$sample");
	my $bin_dir = "$dir_process/$sample/Bin";
	system ("mkdir -m 755 -p $bin_dir") unless (-d $bin_dir);
	system ("mkdir -m 755 -p $dir_result/$sample") unless (-d "$dir_result/$sample");
	system ("mkdir -m 755 -p $dir_sh/$sample") unless (-d "$dir_sh/$sample");	
	my %ass;
	open STAT, "$sample{$sample}{'stat'}";
	<STAT>;
	while (<STAT>)
	{
		chomp;
		my @array = split /\t/;
		$array[0]=~s/Sp\.//;
		$ass{$array[0]} = $array[2];
	}
	close STAT;
	my %stat;
	open CHECKM, "$sample{$sample}{'checkm'}";
	<CHECKM>;
	while(<CHECKM>)
	{
		chomp;
		next if ($_ =~ /^\-\-\-\-/);
		next if ($_ =~ /Bin Id/);
		my @array = split /\s{2,}/;
		$array[1] =~ s/^\s+//g;
		$stat{$array[1]}{'comp'} = $array[-3];
		$stat{$array[1]}{'cont'} = $array[-2];
	}
	my $species_num;
	my @bin_file = `ls $sample{$sample}{'bin'}/*fa`;
	foreach my $bin_file (@bin_file)
	{
		chomp ($bin_file);
		next if (-z $bin_file);
		my $file_name = (split /\//,$bin_file)[-1];
		my $species_id = $file_name;
		$species_id =~ s/\.fa//;	
		next unless (exists $ass{$species_id});
		next unless ($stat{$species_id}{'comp'} > $config{'Comp'} && $stat{$species_id}{'cont'} < $config{'Cont'});
		next if ($ass{$species_id} < $config{'SpeciesLen'});
		$species_num ++;
		open FA, "$bin_file";
		open OUT, ">$bin_dir/Sp.$file_name";
		$/ = "\>";
		<FA>;
		while (<FA>)
		{
			chomp;
			my ($id, $seq) = split /\n+/,$_,2;
			my $ID = (split /\s+/, $id)[0];
			print OUT ">$species_id\_$ID\n$seq";
		}
		$/= "\n";
		close FA;
	}
	print "Sample $sample: $species_num species were used to HGT detecting\n";
	print ASS "$sample\t$species_num\n";
	$shell .= "$Bin/bin/MetaCHIP PI -i $bin_dir -taxon $sample{$sample}{'anno'} -p $sample -r g -t 10 -x fa -force && \\\n";
	$shell .= "$Bin/bin/MetaCHIP BP -p $sample -force -r g $config{'MetaCHIPOpts'} && \\\n";
	$shell .= "tar -zcvf $dir_sh/$sample/$sample\_MetaCHIP_wd/$sample\_g_blastn_results.tar.gz -C $dir_sh/$sample/$sample\_MetaCHIP_wd $sample\_g_blastn_results && \\\n";
	$shell .= "tar -zcvf $dir_sh/$sample/$sample\_MetaCHIP_wd/$sample\_g_blastn_results_filtered_al500bp_cov100.tar.gz -C $dir_sh/$sample/$sample\_MetaCHIP_wd $sample\_g_blastn_results_filtered_al500bp_cov100 && \\\n";
	$shell .= "tar -zcvf $dir_sh/$sample/$sample\_MetaCHIP_wd/$sample\_g_get_SCG_tree_wd.tar.gz -C $dir_sh/$sample/$sample\_MetaCHIP_wd $sample\_g_get_SCG_tree_wd && \\\n";
	$shell .= "mv $dir_sh/$sample/$sample\_MetaCHIP_wd/$sample\_HGT*/$sample\_*_Flanking_region_plots $dir_sh/$sample/$sample\_MetaCHIP_wd/$sample\_g_figures && \\\n";
	$shell .= "tar -zcvf $dir_sh/$sample/$sample\_MetaCHIP_wd/$sample\_g_figures.tar.gz -C $dir_sh/$sample/$sample\_MetaCHIP_wd $sample\_g_figures && \\\n";
	$shell .= "rm -r $dir_sh/$sample/$sample\_MetaCHIP_wd/{$sample\_g_blastn_results,$sample\_g_blastn_results_filtered*,$sample\_g_get_SCG_tree_wd,$sample\_g_blastdb,$sample\_g_figures} && \\\n";
	$shell .= "mv $dir_sh/$sample/$sample\_MetaCHIP_wd $dir_process/$sample/ && \\\n";
	$shell .= "perl $Bin/bin/Blat_for_BM.pl -file $dir_process/$sample/$sample\_MetaCHIP_wd/$sample*/$sample\*_detected_HGTs.txt -bin $bin_dir -out $dir_process/$sample -gene $dir_process/$sample/$sample\_MetaCHIP_wd/$sample\*_prodigal_output/ -id 10 -mis 1000 -gap 1000 && \\\n";
	$shell .= "tar -zcvf $dir_process/$sample/$sample\_MetaCHIP_wd/$sample\_g_prodigal_output.tar.gz -C $dir_process/$sample/$sample\_MetaCHIP_wd $sample\_g_prodigal_output && \\\n";
	$shell .= "rm -r $dir_process/$sample/$sample\_MetaCHIP_wd/$sample\_g_prodigal_output && \\\n";
	$shell .= "tar -zcvf $dir_process/$sample/Bin.tar.gz -C $dir_process/$sample Bin && \\\n";
	$shell .= "tar -zcvf $dir_process/$sample/Blat_result.tar.gz -C $dir_process/$sample Blat_result && \\\n";
	$shell .= "rm -r $dir_process/$sample/{Blat_result,Bin} && \\\n";
	$shell .= "cp $dir_process/$sample/MetaCHIP_Blat_Merge_BM.xls $dir_result/$sample";
	AnaMethod::generateShell("$dir_sh/$sample/MetaCHIP_$sample.sh",$shell);
	$shell = "";
	system("sh $dir_sh/$sample/MetaCHIP_$sample.sh 1>$dir_sh/$sample/MetaCHIP_$sample.sh.o 2>$dir_sh/$sample/MetaCHIP_$sample.sh.e") if ($run =~/y/i && !defined $qopts);
	print DEP "$dir_sh/$sample/MetaCHIP_$sample.sh:$config{'HGTVf'}\n";
	print LIST "$sample\t$dir_process/$sample/$sample\_MetaCHIP_wd/$sample\_HGT*/$sample\_g_detected_HGTs_donor_genes.ffn\t$dir_process/$sample/$sample\_MetaCHIP_wd/$sample\_HGT*/$sample\_g_detected_HGTs_donor_genes.faa\t$dir_process/$sample/MetaCHIP_Blast_Merge.xls\n";	
}
close DEP;
close LIST;
close ASS;

#================ subroutine ===================
sub readconf {
	my ($config) = @_;
	my %conf;
	open IN, "$config" || die "can not open: $config, $!\n";
	while (<IN>)
	{
		chomp;
		next if(/^\s*$/ || /^\s*\#/);
		$_ =~ s/^\s*//;
		$_ =~ s/#(.)*//;
		$_ =~ s/\s*$//;
		if (/^(\w+)\s*=\s*(.*)$/xms) {
			next if ($2 =~ /^\s*$/);
			my $key = $1;
			my $value = $2;
			$value =~ s/\s*$//;
			$conf{$key} = $value;
		}
	}
	return \%conf;
}
#------- THE END ---------
