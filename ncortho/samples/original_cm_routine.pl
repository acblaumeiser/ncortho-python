###################################################################################
#################### SEARCH UP AND DOWN STREAM SEQUENCE OF NCRNA ##################
###################################################################################

# |---DS---|=============|---ncRNA----|=========|---US---| #


	print "Filtering chromosomes and sorting sequences by start position in ROOT GTF file...";
	my @root_keys;
	my %chr_hash;
# 
	foreach(keys %root_gtf_hash){
		my @tmp_array = @{$root_gtf_hash{$_}};  # (47365496 47384273 10 ENST00000585316 ENSG00000265763)
		if 	($tmp_array[2] eq $rna_chr){
			if (not exists $chr_hash{$tmp_array[4]}) {
				$chr_hash{$tmp_array[4]}=\@tmp_array;
			}	
		}
		@root_keys = sort { ${$chr_hash{$a}}[0] <=> ${$chr_hash{$b}}[0]} keys(%chr_hash);
	}
	print "done\n";
	my @US=(); # array of all Upstream genes of ncRNA
	my @DS=(); # array of all Downstream genes of ncRNA
	my @MS=(); # ncRNA matches within Protein ( e.g. intron )
	
	foreach(@root_keys){
		my $ensg = $_;
		my @tmp_array = @{$chr_hash{$ensg}};

		if 	($tmp_array[0] > $rna_stop){
			#print "Stop position of intergenic region $tmp_array[0]\n";
			push (@US,\@tmp_array);
		}
		elsif 	($tmp_array[1] < $rna_start){
			#print "Start position of intergenic region $tmp_array[1]\n";
			push (@DS,\@tmp_array);
		}
		elsif 	(($tmp_array[0] <= $rna_start) && ($tmp_array[1] >= $rna_stop)){
			#print "+++++ INTRON HIT $rna_name : @tmp_array\n";  # debug
			push (@MS,\@tmp_array);
		}
	}
	my $rna_file = $outpath."/core_orthos.fa";
	my $rna_file_seq = "";
	# a test that prints the number of upstream annotations
	#my $test = scalar(@US);
	#print "$test\n";
	if ((scalar(@US) == 0) && (scalar(@MS)== 0)){ 
		print "ERROR : Did not find Upstream Sequence! Terminating\n"; 
		# in the root genome, no gene was found upstream of the miRNA position, therefore a shared syntenic region cannot be constructed
		exit;
	}	
	elsif ((scalar(@DS) == 0) && (scalar(@MS)== 0)){
		print "ERROR : Did not find Downstream Sequence! Terminating\n"; 
		# in the root genome, no gene was found downstream of the miRNA position, therefore a shared syntenic region cannot be constructed
		exit;
	}	

	# routine for the case that the miRNA of the root genome is located within a gene
	elsif (scalar(@MS)> 0){
		my %MS_hash;
		print "ncRNA is located within PROTEIN\n";
		foreach(keys %oma_hash){
			my $species = $_;
			my %tmp_hash = %{$oma_hash{$species}};
			my $ms_index=0;		# if more than one gene is in @MS

			if (exists $tmp_hash{${$MS[$ms_index]}[4]}){
				print "$species has ortholog of ${$MS[$ms_index]}[4]\n";
				$MS_hash{$species} = $tmp_hash{${$MS[$ms_index]}[4]};
				print "==== @{$tmp_hash{${$MS[$ms_index]}[4]}}\n";
			}
			else{
				$ms_index+=1;
				my $found_ms_seq = 0;
				my $ms_size = @MS-1;
				while (($ms_index <= $ms_size) && ($found_ms_seq == 0)){
					if (exists $tmp_hash{${$MS[$ms_index]}[4]}){
						$MS_hash{$species} = $tmp_hash{${$MS[$ms_index]}[4]};
						$found_ms_seq = 1;
					}
					$ms_index+=1;
				}
				if (not exists $MS_hash{$species}){
					print "ncRNA within Protein: No ortholog found for $species\n";
					delete $core_gtf_hash{$species};
                   	delete $core_genome_hash{$species};
				}
			}
		}
		# for each potential coortholog, a blast run is done, then the one with the highest score is selected
		foreach (keys %MS_hash){
			my $species = $_;
			my @tmp_array = @{$MS_hash{$species}};
			print "co ortholog array :@tmp_array\n";
			my %tmp_hash = %{$core_gtf_hash{$species}};
			my $core_genome_file = $core_genome_hash{$species};
			print "$core_genome_file\n";
			my $core_rna_start;
			my $core_rna_stop;
			my $core_rna_chr;
			my $core_rna_strand;
			my $core_rna_score = 0;
			my $highscore_outfile = "";
    	   
			foreach(@tmp_array){
				my $ortholog = $_;
				my @pos = @{$tmp_hash{$ortholog}};
				print "positions: @pos\n"; # 43620095 43652148 1 ENSMMUT00000014574 ENSMMUG00000014970
               	# parse intergenic sequence between up and downstream   
                my $inter_seq   = &genome_parser($core_genome_file,$pos[0],$pos[1],$pos[2],2);
				my $len_seq = length($inter_seq);
               	print "INTERSEQ-LEN: $len_seq\n";
	
	            my $outfile     = $outpath."/".$species."_".$pos[4].".interseq"; # G.gorilla_ENSGGOG00000034972.interseq
	            open(OUTPUT,">",$outfile);
	            print OUTPUT ">$pos[2]\n$inter_seq";
	            close(OUTPUT);
			
				# search for the initial ncRNA within the identified ortholog PCT sequences of each core species
	            system("$formatdb -dbtype nucl -in $outfile"); 
        	   	my $interseq_blast_out = $outfile.".blastout";
				#system("$blastn -p blastn -d $outfile -i $ncRNA -o $interseq_blast_out -m 8"); #perform the blastn search
#                system("$blastn -db $root_genome -query $ncRNA -out $ncPosBLASTout -outfmt 6");
				system("$blastn -task blastn -db $outfile -query $ncRNA -out $interseq_blast_out -outfmt 6 -num_threads $cpu"); #perform the blastn search

	            if (-z $interseq_blast_out){
	            	print "No BLAST HIT found. Skipping\n";
					next;
	            }
	            my ($pred_rna_start,$pred_rna_stop,$pred_rna_chr,$pred_rna_name,$pred_rna_strand,$pred_rna_score) = &blast_parser($interseq_blast_out);

				# incooporate length criterion on the matched region by the blast search
				if (($pred_rna_stop-$pred_rna_start)<($ncRNA_seq_len*$min_seq_len)){
               		print STDERR "§§§§ The potential ortholog of $species was too short.\n";	###!!!### dies here
               	    next;
               	}
				if ($core_rna_score==0){
					$core_rna_start		= $pred_rna_start;
					$core_rna_stop		= $pred_rna_stop;
					$core_rna_chr		= $pred_rna_chr;
					$core_rna_strand	= $pred_rna_strand;
					$core_rna_score		= $pred_rna_score;
					$highscore_outfile 	= $outfile;
				}
				else{
					if($pred_rna_score > $core_rna_score){
				        $core_rna_start         = $pred_rna_start;
                        $core_rna_stop          = $pred_rna_stop;
                        $core_rna_chr           = $pred_rna_chr;
                	    $core_rna_strand        = $pred_rna_strand;
		                $core_rna_score         = $pred_rna_score;
						$highscore_outfile	= $outfile;
					}
				}
			}
			if(not defined $core_rna_score){
		       	print STDERR "No BLASTN HIT found , skipping $species!\n";
			}
			if ($highscore_outfile eq ""){
				next;
			}
            my $pred_rna = &genome_parser($highscore_outfile,$core_rna_start,$core_rna_stop,$core_rna_chr,$core_rna_strand);
			$rna_file_seq.=">$species\n$pred_rna\n";
		}
	}

	# the shared syntenic region does exist in the root genome, now we check if it is present in each of the core species
	# the value in MIP is the number of insertions which are allowed within the syntenic region of the core genome
	else{	
		my %US_hash;
		my %DS_hash;
		print "ncRNA in intergenic region, construct shared syntenic region: US--------miRNA-----------DS\n";
		#my $free_insertions = $max_intra_prot; # allowed insertions of proteins between up and downstream - default 0
		#print "INITIAL NUMBER OF FREE INSERTIONS: $free_insertions\n";
		foreach(keys %oma_hash){ # wird für jede species ausgeführt
			my $us_index = 0;
			my $ds_index = -1;
			my $species = $_;
			print "###################################################################################\n";
			print "######### Construction of shared syntenic region for species:\t$species \n";
			print "###################################################################################\n";
			my %tmp_hash = %{$oma_hash{$species}}; # contains oma orthologs for that species
			my $free_insertions = $max_intra_prot; # allowed insertions of proteins between up and downstream - default 0 ---------------- Mirko edit
			print "INITIAL NUMBER OF FREE INSERTIONS FOR SPECIES $species:\t $free_insertions\n"; #----------------- Mirko edit
		
			print "initial US gene of root species: ${$US[$us_index]}[4] \n";

##################################################################################################
# check the upstream region of the miRNA 
##################################################################################################

			# if the hit is found right away
			if (exists $tmp_hash{${$US[$us_index]}[4]}){
				print "$species: ortholog exists in the first place, no insertions needed\n ";
				$US_hash{$species} = $tmp_hash{${$US[$us_index]}[4]};
			
				print " +++++ print current US hash INITIAL +++++\n";
				my $key = $species;
				print "KEY: $key ## ";
				my @values = @{$US_hash{$species}};
				foreach (@values) {
					print "$_\t";
				}
				print"\n";

			}
		# or else we need to allow insertions
			else{
				$us_index+=1;
				my $found_us_seq = 0;
				my $us_size = @US-1;
				while (($free_insertions > 0) && ($us_index <= $us_size) && ($found_us_seq==0)){
					$free_insertions-=1;
					print"UPDATE ON FREE INSERTIONS WHILE CHECKING US:  $free_insertions\n"; ## testing
					if (exists $tmp_hash{${$US[$us_index]}[4]}){
						$US_hash{$species} = $tmp_hash{${$US[$us_index]}[4]};
						$found_us_seq = 1;
					
						print " +++++ print current US hash UPDATE +++++\n";
						my $key = $species;
						print "KEY: $key ## ";
						my @values = @{$US_hash{$species}};
						foreach (@values) {
						print "$_\t";
						}
						print"\n";

					}
					$us_index+=1;
				}
				# if no insertions are left anymore, cancel routine and delete species from potential core ortholog set, cause shared synteny can not be established.
				if (not exists $US_hash{$species}){
					if ($free_insertions == 0){
						print "UPSTREAM : No free insertions remaining. Did not find orthologous sequence.\n";
					}	
					elsif ( $us_index == scalar(@US)){
						print "UPSTREAM : No potential orthologs left on upstream chromosome\n";
					}
					print "UPSTREAM : No ortholog found for $species! Skipping species!\n";
					delete $core_gtf_hash{$species};
					delete $core_genome_hash{$species};
				}
			}
##################################################################################################
# check the downstream region of the miRNA 
##################################################################################################
		
			print "initial DS gene of root species: ${$DS[$ds_index]}[4] \n";

			# if the hit is found right away
			if (exists $tmp_hash{${$DS[$ds_index]}[4]}){
				$DS_hash{$species} = $tmp_hash{${$DS[$ds_index]}[4]};
				print " +++++ print current DS hash INITIAL +++++\n";
				my $key = $species;
				print "KEY: $key ## ";
				my @values = @{$DS_hash{$species}};
				foreach (@values) {
					print "$_\t";
				}
				print"\n";
			}

			# or else we need to allow insertions		
			else{
				$ds_index-=1;
				my $found_ds_seq = 0;
				my $ds_size = @DS-1;
				while (($free_insertions > 0) && ($ds_index <= $ds_size) && ($found_ds_seq==0)){
					$free_insertions-=1;
					print"UPDATE ON FREE INSERTIONS WHILE CHECKING DS  $free_insertions\n";
					if (exists $tmp_hash{${$DS[$ds_index]}[4]}){
						$DS_hash{$species} = $tmp_hash{${$DS[$ds_index]}[4]};
						$found_ds_seq = 1;
						print " +++++ print current US hash UPDATE +++++\n";
						my $key = $species;
						print "KEY: $key ## ";
						my @values = @{$DS_hash{$species}};
						foreach (@values) {
							print "$_\t";
						}
						print"\n";
					}
					$ds_index+=1;
				}

				# if no insertions are left anymore, cancel routine and delete species from potential core ortholog set, cause shared synteny can not be established.
				if (not exists $DS_hash{$species}){
					if ($free_insertions == 0){
						print "DOWNSTREAM : No free insertions remaining. Did not find orthologous sequence.\n"; 
					}
					elsif ( $ds_index == scalar(@DS)){
						print "DOWNSTREAM : No potential orthologs left on downstream chromosome\n";
					}
					print "DOWNSTREAM : No ortholog found for $species! Skipping species!\n";
					delete $core_gtf_hash{$species};
					delete $core_genome_hash{$species};
				}
			}
		}

		if (!keys %US_hash){
			print STDERR "Absolutely no ortholog UPSTREAM sequences found! Exiting...\n";
			die;
		}
		elsif (!keys %DS_hash){
			print STDERR "Absolutely no ortholog DOWNSTREAM sequences found! Exiting...\n";
			die;
		}
	
	# extract positions (start,stop,chr) from core gtf file based on ENS ids from the constructed shared syntenic regions	

	########################################################################
	#### TO DO : HANDLE CO ORTHOLOGS ( array size > 1 ) ####################
	
		foreach(keys %core_gtf_hash){
			my $species 	= $_;
			my @US_array 	= @{$US_hash{$species}};
			my @DS_array 	= @{$DS_hash{$species}};
		
			my %tmp_hash = %{$core_gtf_hash{$species}};
			if    ((scalar(@US_array) == 1) && (scalar(@DS_array) == 1)){
				$US_hash{$species} = $tmp_hash{$US_array[0]};
				$DS_hash{$species} = $tmp_hash{$DS_array[0]};
			}
			elsif ((scalar(@US_array) == 1) && (scalar(@DS_array)  > 1)){
	        	$US_hash{$species} = $tmp_hash{$US_array[0]};
				my $US_start= ${$tmp_hash{$US_array[0]}}[0];
				my $highscore_stop=0;
				my $highscore_ensg="";
				foreach(@DS_array){
					my $tmp_ensg = $_;
					my $tmp_stop = ${$tmp_hash{$tmp_ensg}}[1];
					if (($highscore_ensg eq "") && ($tmp_stop < $US_start)){
						$highscore_ensg = $tmp_ensg;
						$highscore_stop = $tmp_stop;
					}
					else{
						if (($tmp_stop > $highscore_stop) && ($tmp_stop < $US_start)){
							$highscore_stop = $tmp_stop;
							$highscore_ensg = $tmp_ensg;
						}	
					}
				}
				$DS_hash{$species} = $tmp_hash{$highscore_ensg};
			}
			elsif ((scalar(@US_array)  > 1) && (scalar(@DS_array) == 1)){
	        	$DS_hash{$species} = $tmp_hash{$DS_array[0]};
	            my $DS_stop= ${$tmp_hash{$DS_array[0]}}[0];
	            my $highscore_start=0;
	            my $highscore_ensg="";
	            foreach(@US_array){
	            	my $tmp_ensg = $_;
	            	my $tmp_start = ${$tmp_hash{$tmp_ensg}}[1];
	            	if (($highscore_ensg eq "") && ($tmp_start > $DS_stop)){
	                 	$highscore_ensg = $tmp_ensg;
	                	$highscore_start = $tmp_start;
	            	}
	            	else{
	            		if (($tmp_start > $highscore_start) && ($tmp_start > $DS_stop)){
	            			$highscore_start = $tmp_start;
	            			$highscore_ensg = $tmp_ensg;
	            		}
	            	}
	            }
	            $US_hash{$species} = $tmp_hash{$highscore_ensg};
			}
			elsif ((scalar(@US_array)  > 1) && (scalar(@DS_array) > 1)){
				my %tmp_us_hash; # {ensg}=$start
				my %tmp_ds_hash; # {ensg}=$stop
				foreach (@US_array){
					$tmp_us_hash{${$tmp_hash{$_}}[4]}=${$tmp_hash{$_}}[0];
				}
				foreach (@DS_array){
					$tmp_ds_hash{${$tmp_hash{$_}}[4]}=${$tmp_hash{$_}}[1];
				}
			}
		}

		### check how many proteins are located between up and downstream
		my %count_hash;
		foreach(keys %core_gtf_hash){
			my $species  = $_;
			my $us_start = ${$US_hash{$_}}[0];
			my $ds_stop  = ${$DS_hash{$_}}[1];
			my $ds_chr   = ${$DS_hash{$_}}[2];
			my $counter = 0;
		
			# check if all values needed to construct the blast library are indeed present.
			# a shared syntenic region might be constructed based on the oma (orthology) data, even though a gene used 
			# in that region might not be annotated in the gtf file. In that case the region has to be skipped because of missing data.
	
			if(not defined $us_start){
				print "Error: For species $species Upstream Start Pos. doesn't exist; skipping species cause of missing GTF annotations!\n";
				delete $count_hash{$species};
				delete $US_hash{$species};
				delete $DS_hash{$species};
	            delete $core_genome_hash{$species};
				next;
			}
			elsif (not defined $ds_stop){
				print "Error: For species $species Downstream Stop Pos. doesn't exist; skipping species cause of missing GTF annotations!\n";
				delete $count_hash{$species};
				delete $US_hash{$species};
				delete $DS_hash{$species};
				delete $core_genome_hash{$species};
				next;
			}
			elsif (not defined $ds_chr){
				print "Error: For species $species Chromosome doesn't exist; skipping species cause of missing GTF annotations!\n";
				delete $count_hash{$species};
				delete $US_hash{$species};
				delete $DS_hash{$species};
				delete $core_genome_hash{$species};
				next;
			}
			else { 
				# also added by mirko
				print "species $species\t start $us_start\t stop $ds_stop \t chr $ds_chr\n";
				my %tmp_hash = %{$core_gtf_hash{$species}};
				foreach(keys %tmp_hash){
					my @tmp_array = @{$tmp_hash{$_}};
					if (($tmp_array[0] > $ds_stop) && ($tmp_array[1] < $us_start)&& ($ds_chr eq $tmp_array[2])){
						$counter+=1;
					}
				}
				$count_hash{$species}=$counter;
			}	
		}
		
		# if there is one core species in the core set, that has more inter proteins than the threshold allows (default 0)
		# delete the species and only take the species that fit the criteria
		foreach(keys %count_hash){
			my $species = $_;
			if ($count_hash{$species}>$max_intra_prot){
				print "$species deleted because of intergenic protein coding sequence\n";
				delete $count_hash{$species};
				delete $US_hash{$species};
				delete $DS_hash{$species};
				delete $core_genome_hash{$species};
			}
		}
		# count the species that left after deleting the ones over the threshold
		my @carray=keys(%count_hash);
		my $csize=scalar(@carray); # how many species are left	
		my $take='no';
		if (keys %count_hash){
			foreach(keys %count_hash){
				if ($count_hash{$_} <= $max_intra_prot  ){
					$take='yes';
				}
				else{
					$take='no';
				}
			}
		}
		else{
			print STDERR "All core species over threshold! ";
			die "Terminating!\n";
		}
		### parse the intergenic region in the genome of each core species ###
		foreach(keys %core_genome_hash){	
			my $species	= $_;
			my $genome_file = $core_genome_hash{$species};
			my $start	= ${$DS_hash{$_}}[0];
			my $stop	= ${$US_hash{$_}}[1];
			my $chr		= ${$US_hash{$_}}[2];
			my $strand = 2; # set strand to 2, because for intergen region, it is unimportant
			# parse intergenic sequence between up and downstream
			#print "---> start genome parser with: genome:$genome_file\tstart:$start\tstop:$stop\tchr:$chr\tstrand:$strand\n";

			my $inter_seq	= &genome_parser($genome_file,$start,$stop,$chr,$strand);
			#print "$inter_seq\n";		
			my $len_seq = length($inter_seq);
			print "INTERSEQ-LEN: $len_seq\n";
			my $outfile	= $outpath."/".$species.".interseq";	
			open(OUTPUT,">",$outfile);
			print OUTPUT ">$chr\n$inter_seq";
			close(OUTPUT);
			#system("$formatdb -p F -i $outfile");
			system("$formatdb -dbtype nucl -in $outfile");
			#blastn search in intergenic sequence for ncRNA
			my $interseq_blast_out = $outfile.".blastout";
			system("$blastn -task blastn -db $outfile -query $ncRNA -out $interseq_blast_out -outfmt 6 -num_threads $cpu");

			if (-z $interseq_blast_out){
				print STDERR "No BLASTN HIT found , skipping $species!\n";
				next;
			}
			my ($pred_rna_start,$pred_rna_stop,$pred_rna_chr,$pred_rna_name,$pred_rna_strand,$pred_rna_score) = &blast_parser($interseq_blast_out);
			# if the found sequence has less than a $min_seq_len of the the query RNA, it is discarded
			if (($pred_rna_stop-$pred_rna_start)<($ncRNA_seq_len*$min_seq_len)){
				delete $core_genome_hash{$species};
				print STDERR "§§§§ The potential ortholog of $species was too short\n";
				next;
			}
			my $pred_rna = &genome_parser($outfile,$pred_rna_start,$pred_rna_stop,$pred_rna_chr,$pred_rna_strand);
			$rna_file_seq.=">$species\n$pred_rna\n";
		}

	}	
	if ($rna_file_seq eq ""){
		# if this is inactive, t_coffee might run into a core dump, besides... alignment with  1 seq does not make sense...
		print STDERR "§§ The sequence file for the alignment only contains 1 sequence - Exiting $ncRNA!\n";
		die;
	}

	# generates core orthologs outfile 
	$rna_file_seq.=">query\n$ncRNA_seq";  
	open(RNAOUT,">",$rna_file);
	print RNAOUT $rna_file_seq;
	close(RNAOUT);

	print "STEP 01\n";
	system("$tcoffee -cpu $cpu -special_mode=rcoffee -in $rna_file -output=clustalw_aln > $alignment"); 

	print "STEP 02\n";
	system("$tcoffee -other_pg seq_reformat -in $alignment -action +add_alifold -output stockholm_aln > $stockholm_aln");

	print "STEP 03\n";
	system("$cmbuild -F $covariance_model $stockholm_aln");
#	system("$cmbuild $covariance_model $stockholm_aln");
	print "STEP 04\n";
	system("$cmcalibrate --cpu $cpu $covariance_model");
}
