package AnaMethod;

my $package_path = "/share/backup/qiukl/pipeline/packages";

########################################################################
#Abstract:
#	This Function is use to generate script.
#para:
#	$output_shell: Script to be output
#	$content: Main content in the script
########################################################################
sub generateShell
{
	my ($output_shell, $content, $finish_string) = @_;
	unlink glob "$output_shell.*";
	$finish_string ||= "Still_waters_run_deep";
	open OUT,">$output_shell" or die "Cannot open file $output_shell:$!";
	print OUT "#!/bin/bash\n";
	print OUT "echo ==========start at : `date` ==========\n";
	print OUT "$content && \\\n";
	print OUT "echo ==========end at : `date` ========== && \\\n";
	print OUT "echo $finish_string 1>&2 && \\\n";
	print OUT "echo $finish_string > $output_shell.sign\n";
	close OUT;

}
#########################################################################
#Abstract:
#  This function is use to generate scripts of bwa aln
#para:
#	$output_shell: Script to be output
#	$bwa: path of bwa
#	$bwa_aln_para: parameters of  bwa aln
#	$sai_file_out: sai file to be output
#	$input_fq_file : input fastq file 
#########################################################################
sub bwaAln
{
	my ($bwa, $bwa_aln_para, $sai_file_out, $ref, $input_fq_file) = @_;
	my $content = "$bwa aln $bwa_aln_para -f $sai_file_out $ref $input_fq_file";
	return $content;
}

########################################################################################
#Abstract:
#  This function is use to generate scripts of bwa sampe and samtools sort
#para:
########################################################################################
sub bwaSampeAndSort
{
	my ($bwa, $samtools, $bwa_sampe_para, $ref, $sai1, $sai2, $fq1, $fq2, $output_prefix) = @_;
	my $content = "$bwa sampe $bwa_sampe_para $ref $sai1 $sai2 $fq1 $fq2 | $samtools view -Sb -T $ref - >$output_prefix.temp.bam && \\\n";
	$content .= "$samtools sort -m 3500000000 $output_prefix.temp.bam $output_prefix && \\\n";
	$content .= "rm -f $output_prefix.temp.bam && \\\n";
	$content .= "$samtools index $output_prefix.bam";
	
	return $content;
}

sub soap3dp
{
	my ($soap3dp, $novosort, $soap3dp_para, $novosort_para, $ref, $fq1, $fq2, $lane, $lib, $pu, $output_prefix,$sortTmpdir) =@_;
	my $content = "$soap3dp pair $ref $fq1 $fq2 -D $lane -R \'LB:$lib\\tPL:Illumina\\tPU:$pu\' $soap3dp_para -o $output_prefix && \\\n";
	$content .= "$novosort $novosort_para -t $sortTmpdir -o $output_prefix.bam ";
	for (my $i=1; $i<=10;$i++) {
		$content .= "$output_prefix.gout.$i ";
	}
	$content .= "$output_prefix.dpout.1 $output_prefix.unpair";

	return $content;
}

sub bwaMem
{
	my ($bwa, $bwa_mem_para, $samtools, $view_para, $sort_para, $ref, $fq1, $fq2, $output_prefix) = @_;
	my $content = "$bwa $bwa_mem_para $ref $fq1 $fq2 |$samtools view $view_para - -o $output_prefix.temp.bam && \\\n";
	$content .= "$samtools sort $sort_para $output_prefix.temp.bam $output_prefix";

	return $content;
}

########################################################################################
#Abstract:
#  This function is use to generate scripts of bwa sampe and picard/SortSam.jar sort 
#  and samtools index
#para:
########################################################################################
sub bwaSampePicardSortSamtoolsIndex
{
	my ($java, $bwa, $picard_jar, $samtools, $bwa_sampe_para, $ref, $sai1, $sai2, $fq1, $fq2, $java_tmp, $tmp, $output_prefix, $simple) = @_;
	my $content = "$bwa sampe $bwa_sampe_para $ref $sai1 $sai2 $fq1 $fq2 | $samtools view -Sb -T $ref - >$output_prefix.temp.bam && \\\n";
	if($simple == 1)
	{
		$content .= "$java -Xmx3g -Djava.io.tmpdir=$java_tmp -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $picard_jar/SortSam.jar I=$output_prefix.temp.bam O=$output_prefix.bam TMP_DIR=$tmp SO=coordinate VALIDATION_STRINGENCY=SILENT && \\\n";
	}
	else
	{
		$content .= "$java -Xmx3g -Djava.io.tmpdir=$java_tmp -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $picard_jar/FixMateInformation.jar I=$output_prefix.temp.bam O=$output_prefix.fixmate.bam TMP_DIR=$tmp SO=coordinate VALIDATION_STRINGENCY=SILENT && \\\n";
		$content .= "mv $output_prefix.fixmate.bam $output_prefix.temp.bam && \\\n";
		$content .= "$samtools calmd -b $output_prefix.temp.bam $ref > $output_prefix.bam && \\\n";
	}
	$content .= "rm -f $output_prefix.temp.bam && \\\n";
	$content .= "$samtools index $output_prefix.bam";
	
	return $content;
}


########################################################################################
#Abstract:
#	This function is use to generate scripts of merge bams by picard
#para:
#	$output_shell: Script to be output
#	$merge_sam_tool: the tool use to merge bam files
#	$bam_array: reference of sorted bam array
#	$output_bam: merged output bam
#	$tmp_dir: temp directory
#	$rm_input_bams: 1 or 0, 1 means remove input bam files
#
########################################################################################
sub mergeBamsByPicard
{
	my ($output_shell, $java, $merge_sam_tool, $bam_array, $output_bam, $java_tmp, $tmp_dir, $rm_input_bams) = @_;
	my $lane_number = scalar @{$bam_array};	

	my $content;
	my $bam_str;
	if ($lane_number > 1) {
		$content = "$java -Xmx2G -Djava.io.tmpdir=$java_tmp -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $merge_sam_tool ";
		foreach my $inbam (@{$bam_array}) {
			$content .= "I=$inbam ";
			$bam_str .= "$inbam $inbam.bai ";
		}
		$content .= "O=$output_bam TMP_DIR=$tmp_dir SO=coordinate AS=true VALIDATION_STRINGENCY=SILENT";
	}
	else {
		$content = "cp ${$bam_array}[0] $output_bam";
		$bam_str = "${$bam_array}[0] ${$bam_array}[0].bai ";
	}

	if ($rm_input_bams) {
		$content .= " && \\\n rm -f $bam_str";
	}

	generateShell($output_shell, $content);
}


########################################################################################
#Abstract:
#	This function is use to generate scripts of merge bams by picard and samtools index
#para:
#	$output_shell: Script to be output
#	$merge_sam_tool: the tool use to merge bam files
#	$samtools: the tool use to build index for the bam file
#	$bam_array: reference of sorted bam array
#	$output_bam: merged output bam
#	$tmp_dir: temp directory
#	$rm_input_bams: 1 or 0, 1 means remove input bam files
#
########################################################################################
sub mergeBamsByPicardSamtoolsIndex
{
	my ($output_shell, $java, $merge_sam_tool, $samtools, $bam_array, $output_bam, $java_tmp, $tmp_dir, $rm_input_bams) = @_;
	my $lane_number = scalar @{$bam_array};	

	my $content;
	my $bam_str;
	if ($lane_number > 1) {
		$content = "$java -Xmx2G -Djava.io.tmpdir=$java_tmp -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $merge_sam_tool ";
		foreach my $inbam (@{$bam_array}) {
			$content .= "I=$inbam ";
			$bam_str .= "$inbam $inbam.bai ";
		}
		$content .= "O=$output_bam TMP_DIR=$tmp_dir SO=coordinate AS=true VALIDATION_STRINGENCY=SILENT";
	}
	else {
		$content = "cp ${$bam_array}[0] $output_bam";
		$bam_str = "${$bam_array}[0] ${$bam_array}[0].bai ";
	}

	if ($rm_input_bams) {
		$content .= " && \\\n rm -f $bam_str";
	}

	$content .= " && \\\n $samtools index $output_bam ";


	generateShell($output_shell, $content);
}


########################################################################################
#Abstract:
#	This function is use to mark PCR duplicates by picard
#para:
#	$output_shell: Script to be output
#	$rm_input_bams: 1 or 0, 1 means remove input bam files
#
#######################################################################################
sub markDuplicatesByPicard
{
	my ($output_shell, $java, $markdup_tool, $inbam, $outbam, $java_tmp, $tmp_dir, $rm_input_bams) = @_;

	my $content = "$java -Xmx2G -Djava.io.tmpdir=$java_tmp -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $markdup_tool REMOVE_DUPLICATES=false I=$inbam O=$outbam METRICS_FILE=$outbam.mat TMP_DIR=$tmp_dir VALIDATION_STRINGENCY=SILENT ";
	if ($rm_input_bams) {
		$content .= " && \\\n rm -f $inbam $inbam.bai";
	}
	generateShell($output_shell, $content);
}

########################################################################################
#Abstract:
#	This function is use to mark PCR duplicates by picard ,then build index for the new bam by samtools
#para:
#	$output_shell: Script to be output
#	$samtools: path of samtools
#	$rm_input_bams: 1 or 0, 1 means remove input bam files
#
#######################################################################################
sub markDuplicatesByPicardSamtoolsIndex
{
	my ($output_shell, $java, $markdup_tool, $samtools, $inbam, $outbam, $java_tmp, $tmp_dir, $rm_input_bams) = @_;

	my $content = "$java -Xmx2G -Djava.io.tmpdir=$java_tmp -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar $markdup_tool REMOVE_DUPLICATES=false I=$inbam O=$outbam METRICS_FILE=$outbam.mat TMP_DIR=$tmp_dir VALIDATION_STRINGENCY=SILENT && \\\n";
	$content .= "$samtools index $outbam";
	if ($rm_input_bams) {
		$content .= " && \\\n rm -f $inbam $inbam.bai";
	}
	generateShell($output_shell, $content);
}

########################################################################################
#Abstract:
#	This function is use to generate scripts of GATK realignment
#para:
#	$output_shell: Script to be output
#	$GATK_tool: path of GATK
#	$ref: reference that GATKrealigner align to
#	$intervals:  One or more genomic intervals over which to operate
#	$target_intervals:
#	$known_sites_array: array of known indel snp sites file
#	$$inbam: input sorted,duplication marked bam
#	$outbam: output realigned bam
#	$java_tmp: template directory for java
#	$rm_input_bams: 1 or 0, 1 means remove input bam files
#
########################################################################################
sub realignByGATK
{
	my ($output_shell, $GATK_tool, $ref, $intervals, $target_intervals, $known_sites_array, $inbam, $outbam, $java_tmp, $rm_input_bams) = @_;
	my ($known_sites, $ks);
	foreach my $ks (@{$known_sites_array}) {
		$known_sites .= "-known $ks ";
	}
	
	my $content = "java -Xmx1g -Djava.io.tmpdir=$java_tmp -jar $GATK_tool -T RealignerTargetCreator -l INFO -I $inbam -R $ref -L $intervals -o $target_intervals $known_sites && \\\n";
	$content .= "java -Xmx4g -Djava.io.tmpdir=$java_tmp -jar $GATK_tool -T IndelRealigner -l INFO -I $inbam -R $ref -L $intervals -targetIntervals $target_intervals $known_sites -o $outbam"; 

	if ($outbam =~ /(\S+)\.bam$/) {
		my $prefix = $1;
		$content .= " && \\\n mv $prefix.bai $outbam.bai";
	}
	if ($rm_input_bams) {
		$content .= " && \\\n rm -f $inbam $inbam.bai";
	}

	generateShell($output_shell, $content);
	
}

########################################################################################
#Abstract:
#	This function is use to generate scripts of GATK base recalibration
#para:
#	$output_shell: Script to be output
#	$GATK_tool: path of GATK
#	$ref: reference that GATK align to
#	$intervals:
#	$known_sites_array: array of known indel snp sites file
#	$inbam: input sorted,duplication marked,GATK aligned bam
#	$outbam: output base recalibrated bam
#	$java_tmp: template directory for java
#	$rm_input_bams: 1 or 0, 1 means remove input bam files
#
########################################################################################
sub baseRecalibrateByGATK
{
	my ($output_shell, $GATK_tool, $ref, $intervals, $known_sites_array, $inbam, $outbam, $java_tmp, $rm_input_bams) = @_;
	my ($known_sites, $ks);
	foreach my $ks (@{$known_sites_array}) {
		$known_sites .= "-knownSites $ks ";
	}
	my $content = "java -Xmx4g -Djava.io.tmpdir=$java_tmp -jar $GATK_tool -T BaseRecalibrator -l INFO -I $inbam -R $ref $known_sites -o $inbam.recal_data.table -L $intervals -nct 8 && \\\n";
	$content .= "java -Xmx4g -Djava.io.tmpdir=$java_tmp -jar $GATK_tool -T PrintReads -l INFO -I $inbam -R $ref -BQSR $inbam.recal_data.table -o $outbam -L $intervals"; 

	if ($outbam =~ /(\S+)\.bam$/) {
		my $prefix = $1;
		$content .= " && \\\n mv $prefix.bai $outbam.bai";
	}
	if ($rm_input_bams) {
		$content .= " && \\\n rm -f $inbam $inbam.bai";
	}

	generateShell($output_shell, $content);
}

########################################################################################
#Abstract:
#	This function is use to generate scripts of GATK variation detect.
#para:
#	$output_shell: Script to be output
#	$GATK_tool: path of GATK
#	$ref: reference that GATK align to
#	$intervals:
#	$number_of_threads:
#	$dbsnp:  dbsnp file
#	$inbam: input sorted,duplication marked,GATK aligned bam
#	$outvcf: output variation file 
#	$java_tmp: template directory for java
#
########################################################################################
sub callSmallVariationByGATK
{
	my ($output_shell, $GATK_tool, $ref, $intervals, $number_of_threads, $dbsnp, $inbam, $outvcf, $java_tmp) = @_;
	my $content = "java -Xmx4g -Djava.io.tmpdir=$java_tmp -jar $GATK_tool -T UnifiedGenotyper -l INFO -I $inbam --dbsnp $dbsnp -stand_call_conf 30 -stand_emit_conf     10 -dcov 1000 -L $intervals -R $ref -o $outvcf -nt $number_of_threads -glm BOTH";
	generateShell($output_shell, $content);
}

########################################################################################
#Abstract:
#	This function is use to generate scripts of GATK variation detect.
#para:
#	$output_shell: Script to be output
#	$GATK_tool: path of GATK
#	$ref: reference that GATK align to
#	$hapmap_vcf:  hapmap variation file
#	$_1000G_omni:  1000G omni file
#	$dbsnp:  dbsnp file
#	$invcf: input raw variation file 
#	$outvcf: output variation file 
#	$java_tmp: template directory for java
#
########################################################################################
sub variantRecalibrateByGATK
{
	my ($output_shell, $GATK_tool, $ref, $hapmap_vcf, $_1000G_omni, $dbsnp, $invcf, $outvcf, $java_tmp) = @_;
	my $content = "java -Xmx4g -Djava.io.tmpdir=$java_tmp -jar $GATK_tool -T VariantRecalibrator -l INFO -R $ref -input $invcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_vcf -resource:omni,known=false,training=true,truth=false,prior=12.0 $_1000G_omni -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $dbsnp -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ -recalFile $outvcf.recal -tranchesFile $outvcf.tranches -rscriptFile $outvcf.plot.R --maxGaussians 4 && \\\n";
	$content .= "java -Xmx3g -Djava.io.tmpdir=$java_tmp -jar $GATK_tool -T ApplyRecalibration -l INFO -R $ref -input $invcf --ts_filter_level 99.0 -tranchesFile $outvcf.tranches -recalFile $outvcf.recal -o $outvcf";
	generateShell($output_shell, $content);
}


########################################################################################
#Abstract:
#	This function is use to generate scripts of variation annotation.
#para:
#	$output_shell: Script to be output
#	$vcf : snp or indel vcf file
#	$buildver: --buildver of summarize_annovar.pl, hg18 or hg19
#	$verdbsnp: --verdbsnp of summarize_annovar.pl, 129,130, 132, ..., 137
#	$snp_or_indel: variation type: snp or indel
#	$sample: sample name
#	$output_prefix: output prefix 
#
########################################################################################
sub annotateVariationsByAnnovar
{
	my ($output_shell, $vcf, $buildver, $verdbsnp, $snp_or_indel, $sample, $output_prefix) = @_;
	my $annovar_dir = "$package_path/bin/annovar";
	my $Annovar_stat_for_all = "$package_path/bin/Annovar_stat_for_all.pl";
	my $content = "export PATH=\$PATH:$annovar_dir && \\\n";
	$content .= "perl $annovar_dir/convert2annovar.pl --includeinfo --format vcf4 $vcf > $output_prefix && \\\n";
	$content .= "perl $annovar_dir/summarize_annovar.pl --remove --buildver $buildver --verdbsnp $verdbsnp $output_prefix $annovar_dir/humandb_$buildver && \\\n";
	$content .= "perl $Annovar_stat_for_all -i $output_prefix.genome_summary.csv -o $output_prefix -s $sample -v $snp_or_indel";
	generateShell($output_shell, $content);
}


########################################################################################
#Abstract:
#	This function is use to generate scripts of qsub
#para:
#	$output_shell: Script to be output
#	$task_monitor: path of task monitor program
#	$q: -q of qsub
#	$P: -P of qsub
#	$config: file of inter-dependent relationships(output)
#	$packages: path of python packages
#
########################################################################################
sub outputQsubScript
{
	my ($output_shell, $task_monitor, $q, $P, $config, $packages) = @_;
	open QSUB,">","$output_shell" or die "Cannot open file $output_shell:$!\n";
	print QSUB "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/opt/blc/python-2.6.5/lib\n";
	print QSUB "export PYTHONPATH=\$PYTHONPATH:$packages\n";
	print QSUB "$task_monitor -q $q -i $config -P $P -t 20\n";
	close QSUB;
}
########################################################################################
#Abstract:
#	This function is use to get chromosomes from the fai file of reference.
#para:
#	ref: path of reference
########################################################################################
sub getChromosomesFromRefFai
{
	my ($ref) = @_;
	my @chr;
	open FAI,"$ref.fai" or die "Cannot read the fai file $ref.fai: $!";
	while(<FAI>) {
		chomp;
		my $ccc=(split /\t/,$_)[0];
		next if($ccc =~ /_gl/ || $ccc =~ /_hap/ || $ccc =~ /\./ || $ccc =~ /_random/ || $ccc =~ /chrUn/);
		push @chr,$ccc;
	}
	close FAI;
	return \@chr;
}

########################################################################################
#Abstract:
#	This function is use to get chromosomes from the bam file
#para:
#	$samtools: path of samtools
#	$bam: path of the bam file
########################################################################################
sub getChromosomesFromBam
{
	my ($samtools, $bam) = @_;
	my @chr;
	open BAM,"$samtools view -H $bam |" or die "Cannot read $bam: $!";
	while(<BAM>) {
		chomp;
		next unless ($_ =~ /^\@SQ/);
		my @line = split;
		$line[1] =~ /SN:(\w+)/;
		my $ccc=$1;
		next if($ccc =~ /_gl/ || $ccc =~ /_hap/ || $ccc =~ /\./ || $ccc =~ /_random/ || $ccc =~ /chrUn/);
		push @chr,$ccc;
	}
	close BAM;
	return \@chr;
}

########################################################################################
#Abstract:
#	This function is use to generate scripts of variation annotation.
#para:
#	$annovar_dir: annovar directory
#	$buildver: --buildver of summarize_annovar.pl, hg18 or hg19
#	$verdbsnp: --verdbsnp of summarize_annovar.pl, 129,130, 132, ..., 137
#	$output_prefix: output prefix(format *.annovar) 
#
########################################################################################
sub annotationByAnnovar
{
	my ($annovar_dir, $buildver, $verdbsnp, $output_prefix) = @_;
	my $content = "export PATH=\$PATH:$annovar_dir && \\\n";
	$content .= "perl $annovar_dir/summarize_annovar.pl --remove --buildver $buildver ";
	$content .= "--verdbsnp $verdbsnp " if($buildver =~ /hg19/i);
	$content .= "$output_prefix $annovar_dir/humandb_$buildver";
	return $content;
}


########################################################################################
#Abstract:
#	This function is use to generate scripts of variation annotation for mm9
#para:
#	$annovar_mm9: grogram to annotate mm9
#	$annovar_dir: annovar directory
#	$output_prefix: output prefix(format *.annovar) 
#	$mm9db: annotataion database of mm9
#
########################################################################################
sub annotationByAnnovarForMm9
{
	my ($annovar_mm9, $annovar_dir, $mm9db, $output_prefix) = @_;
	my $content = "export PATH=\$PATH:$annovar_dir && \\\n";
	$content .= "perl $annovar_mm9 -remove -BIN $annovar_dir $output_prefix $mm9db";
	return $content;

}


#########################################################################################
#Abstract:
#	
#Para:
#	$ganno2_db:Optional database. Separated by comma. Available options are [HGMD COSMIC ENCODE YH ESP], default is [ESP,HGMD,COSMIC,YH,ENCODE].  If none are needed, use ' '.
#	$buildver: --buildver of summarize_annovar.pl, hg18 or hg19
#	$annovar_db_new: set the version of dbSNP for annodb.pl (-ganno2 mode), [v138_hg19 v137_hg19]	
#	$type: variation type:snp,indel,snv
#	$output_prefix: output prefix(format *.annovar) 
#
#########################################################################################
sub annotationByAnnodb
{
	my ($annodb, $ganno2_db, $type, $annovar_db_new, $buildver, $output_prefix) = @_; 
	my $content;
	if($type =~ /snp/ || $type =~ /snv/ || $type =~ /indel/) {
		if($type =~ /snp/ || $type =~ /snv/) {
			$content = "perl $annodb --mode snp --if annovar --optdb \'$ganno2_db\' --remove --verdbsnp $annovar_db_new ";
		}
		else {
			$content = "perl $annodb --mode indel --if annovar --optdb \'$ganno2_db\' --remove --verdbsnp $annovar_db_new ";
		}
		if($anno =~ /mm9/i)	{
			$content .= "--buildver $buildver --outfile $output_prefix $output_prefix ";
		}
		else {
			$content .= "--buildver UCSC$buildver --outfile $output_prefix $output_prefix ";
		}
	}
	return $content;

}
##########################################################################################
#Abstract:
#	call small variations(snp/indel) by samtools mpileup for whole exome sequencing data
#parameters:
#	$output_shell: [sh file] output shell
#	$samtools: [tool] path of samtools
#	$sam_mpile_para: [para] parameter string of samtools mpileup
#	$bcftools: [tool] path of bcftools
#	$vcfutils: [tool] path of vcfutils
#	$vcfutils_para: [para] parameter string of vcfutils
#	$filt_vcf: [tool] path of filt_vcf.pl
#	$filt_vcf_para: [para] parameter  of filt_vcf.pl
#	$ref: [file] path of reference
#	$bed: [file] path of target region file(bed format)
#	$bam: [file] alignment result file(bam format)
#	$output_vcf: [file] output variation file (vcf format)
#	$output_filt_vcf: [file] output filtered variation file [vcf format]
##########################################################################################
sub callSnpIndelForWESByMpileup
{
	my ($output_shell, $samtools, $sam_mpile_para, $bcftools, $vcfutils, $vcfutils_para, $filt_vcf, $filt_vcf_para, $ref, $bed, $bam, $output_bcf, $output_vcf, $output_filt_vcf, $chr) = @_; 
	my $content = "$samtools mpileup $sam_mpile_para ";
	if ($chr) {
		$content .= "-l $bed -r $chr -f $ref $bam \\\n";
	}
	else {
		$content .= "-l $bed -f $ref $bam \\\n";
	}
	$content .= "| $bcftools view -bvcg - > $output_bcf && \\\n";
	$content .= "$bcftools view $output_bcf \\\n";
	$content .= "| perl $vcfutils varFilter $vcfutils_para >$output_vcf && \\\n";
	$content .= "perl $filt_vcf $filt_vcf_para -i $output_vcf -o $variation_outdir1/$output_filt_vcf";
	generateShell($output_shell, $content);
}

##############################################################################################
#Abstract:
#	call small variations(snp/indel) by samtools mpileup for whole genome sequencing data
#parameters:
#	$output_shell: [sh file] output shell
#	$samtools: [tool] path of samtools
#	$sam_mpile_para: [para] parameter string of samtools mpileup
#	$bcftools: [tool] path of bcftools
#	$vcfutils: [tool] path of vcfutils
#	$vcfutils_para: [para] parameter string of vcfutils
#	$filt_vcf: [tool] path of filt_vcf.pl
#	$filt_vcf_para: [para] parameter  of filt_vcf.pl
#	$ref: [file] path of reference
#	$bam: [file] alignment result file(bam format)
#	$output_vcf: [file] output variation file (vcf format)
#	$output_filt_vcf: [file] output filtered variation file [vcf format]
#	$chr: [string] set the chromosome name if you only want to call one chromosome's snp/indel
###############################################################################################
sub callSnpIndelForWGSByMpileup
{
	my ($output_shell, $samtools, $sam_mpile_para, $bcftools, $vcfutils, $vcfutils_para, $filt_vcf, $filt_vcf_para, $ref, $bam, $output_bcf, $output_vcf, $output_filt_vcf, $chr) = @_;
	my $content = "$samtools mpileup $sam_mpile_para -f $ref ";
	#$content .= "-l $bed -f $ref $bam \\\n";
	if ($chr) {
		$content .= "-r $chr ";
	}
	$content .= "$bam ";
	$content .= "| $bcftools view -bvcg - > $output_bcf && \\\n";
	$content .= "$bcftools view $output_bcf \\\n";
	$content .= "| perl $vcfutils varFilter $vcfutils_para >$output_vcf && \\\n";
	$content .= "perl $filt_vcf $filt_vcf_para -i $output_vcf -o $output_filt_vcf";
	generateShell($output_shell, $content);
}


sub callSomaticSmallVariationForWGSByVarScan 
{
	my ($output_shell, $samtools, $varscan, $var_para, $ref, $min_baseQ, $normal_bam, $tumor_bam, $java_tmp, $output_snp, $output_indel, $finish_string, $chr) = @_; 
	`rm -f $output_shell*` if(-e "$output_shell");
	$finish_string ||= "Still_waters_run_deep";
	open OUT,">$output_shell" or die $!;
	print OUT "#!/bin/bash\necho ==========start at : `date` ==========\n";
	print OUT "java -Xmx2G -Djava.io.tmpdir=$java_tmp -jar $varscan somatic ";
	print OUT "<($samtools mpileup -f $ref -Q $min_baseQ ";
	print OUT "-r $chr " if (defined $chr);
	print OUT "$normal_bam) <($samtools mpileup -f $ref -Q $min_baseQ ";
	print OUT "-r $chr " if (defined $chr);
	print OUT "$tumor_bam) ";
	if ($min_baseQ == 0) {
		print OUT "$var_para --min-avg-qual 0 --output-snp $output_snp --output-indel $output_indel && \\\n";
	}
	else {
		print OUT "$var_para --output-snp $output_snp --output-indel $output_indel && \\\n";
	}
	print OUT "line=`wc -l $output_snp | awk '{print \$1}'` && \\\n";
	print OUT "echo ==========end at : `date` ========== && \\\n";
	print OUT "if [ \$line -lt 2 && \"$chr\" != \"chrM\"];then echo invalid!;else echo $finish_string 1>&2;";
	print OUT "echo $finish_string > $output_shell.sign;fi\n";
	close OUT;
}

sub callSomaticSmallVariationForWESByVarScan 
{
	my ($output_shell, $samtools, $varscan, $var_para, $ref, $min_baseQ, $normal_bam, $tumor_bam, $java_tmp, $output_snp, $output_indel, $region, $finish_string, $chr) = @_; 
	`rm -f $output_shell*` if(-e "$output_shell");
	$finish_string ||= "Still_waters_run_deep";
	open OUT,">$output_shell" or die $!;
	print OUT "#!/bin/bash\necho ==========start at : `date` ==========\n";
	print OUT "java -Xmx2G -Djava.io.tmpdir=$java_tmp -jar $varscan somatic ";
	print OUT "<($samtools mpileup -f $ref -Q $min_baseQ ";
	print OUT "-l $region ";
	print OUT "-r $chr " if (defined $chr);
	print OUT "$normal_bam) <($samtools mpileup -f $ref -Q $min_baseQ ";
	print OUT "-l $region ";
	print OUT "-r $chr " if (defined $chr);
	print OUT "$tumor_bam) ";
	if ($min_baseQ == 0) {
		print OUT "$var_para --min-avg-qual 0 --output-snp $output_snp --output-indel $output_indel && \\\n";
	}
	else {
		print OUT "$var_para --output-snp $output_snp --output-indel $output_indel && \\\n";
	}
	print OUT "line=`wc -l $output_snp | awk '{print \$1}'` && \\\n";
	print OUT "echo ==========end at : `date` ========== && \\\n";
	print OUT "if [ \$line -lt 2 ];then echo invalid!;else echo $finish_string 1>&2;";
	print OUT "echo $finish_string > $output_shell.sign;fi\n";
	close OUT;
}

sub readConfigureFile
{
	my ($file, $config_hash) = @_;
	open FIN, $file or die "Cannot open file $file:$!\n";
	while (<FIN>)
	{
		chomp;
		next if(/^\s*$/ || /^\s*\#/);
		$_ =~ s/^\s*//;
		$_ =~ s/\s*$//;
		if (/^(\w+)\s*=\s*(.*)$/xms) {
			next if ($2 =~ /^\s*$/);
			my $key = $1;
			my $value = $2;
			$value =~ s/\s*$//;
			${$config_hash}{$key} = $value;
		}
	}
}
	
1;
