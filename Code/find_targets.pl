#########################################
# This file is part of targetFinder.
#
# targetFinder is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# targetFinder is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with targetFinder.  If not, see <https://www.gnu.org/licenses/>.
#
#########################################
# Subroutine findTargets: This script will find potential target locations for 
#  insertion sites. The script takes as input the target genome of interest and
#  outputs a list of suitable locations. The criteria for a suitable 
#  target locations are:
#    - near center of 500bp chunk between 2 convergently transcribed genes (either
#	annotated or predicted as a potential ORF)
#    - on the chromosome (not plasmid)
#    - no repetitive structure within intergenic space
#    - no large repetitive elements (>200bp) within 10kb
# Inputs:
#    - genome sequence to find targets for
# Outputs:
#    - array of sequence positions for the ranked list of targets
#########################################
sub findTargets {
use Bio::Tools::Run::StandAloneBlastPlus;
#########################################
# Interpret input vars
#########################################
#General Parameters
my $genomeFile=@_[0];		#name of the file with info on the current genome
my $v=@_[1]; 			#verbose output off or on
my $log=@_[2];			#write to log file on or off
#Target Parameters
my $minGap=@_[3];		#minimum size of gap between convergent CDS's
my $maxGap=@_[4];		#maximum size of gap between convergent CDS's
my $homRecombSize=@_[5];	#size of flanking regions on either side of 
				# target insertion point 
my $locRep_Nmer=@_[6];		#minimum size of repeated nmers to qualify
				# as a local repitive region (see locRep_Reps
my $locRep_Reps=@_[7];		#minimum number of repeated nmers to qualify
				# as a local repitive region (see locRep_Reps)
my $ORFbuffer=@_[8];		#the size outside of the gap that should be 
				# scanned for potential ORFs (note that this is
				# in addition to scanning the entire gap region,
				# which can be of a size defined above)
my $minORFsize=@_[9];		#minimum size of a potential ORF to be considered
				# a problem
my $maxORFsize=@_[10];		#maximum size of a potential ORF to be considered
				# a problem
my $lgRepBound=@_[11];		#range to search within for large repetitive
				# elements
my $lgRepMinSize=@_[12];	#minimum size of large repetitive elements
my $lgRepMaxSize=@_[13];	#maximum size of large repetitive elements

my $minGapIdentityWarnLen=@_[14]; #minimum number of bases for gap to genome
                                        # identity to generate a warning

my $minIdentityWarn=@_[15];	        #minimum percent identity needed to generate a 
                                        #warning

my $pamSiteSeq=@_[16];	        	#sequence used to identify a PAM site

my $guideRNAlen=@_[17];	                #length	of the guide RNA sequence

my $maxGuideRNAmatches=@_[18];       #max number of base matches to accept a guide RNA

my $smallRNAseqWin=@_[19];	        #small window size to calculate	depth of coverage
my $largeRNAseqWin=@_[20];	        #large window size to calculate depth of coverage

#open genome file
#update log and user on progress
update("Finding targets for genome file: ".$genomeFile."\n",$v,$log);
my $genome=Bio::SeqIO->new(-id=>'curGenome',-file=>'targetGenomes/'.$genomeFile);
my $seq=$genome->next_seq;

#########################################
#Pull coding sequences from genome file
#########################################
#update log and user on progress
update("\tPulling CDS features from genome file...",$v,$log);
#declare variables
my @starts;
my @ends;
my @strands;
my @otherFeatStarts;
my @otherFeatEnds;
#go through each feature and store starts/stops for cds's on each strand as
# appropriate
foreach my $feat ($seq->all_SeqFeatures) {
	#grab the start/end/strand of the feature
	my $curStart=$feat->start;
	my $curEnd=$feat->end;
	my $curStrand=$feat->strand;
	#check if feature is a coding sequence / open reading frame
	if ($feat->primary_tag =~ /cds/i) {
#	if ($feat->primary_tag =~ /gene|cds|orf/i) {
		#-------------------------
		# Scan to make sure it's not a duplicate entry (e.g. there are
		# separate annotations for a CDS and gene)
		#-------------------------
		#if this is the first CDS found, add it
		if (@starts==0) {
			@starts[0]=$curStart;
			@ends[0]=$curEnd;
			@strands[0]=$curStrand;
		}
		my $i=0;
		my $flag=0;
		foreach my $entry (@starts) {
			#if it is a duplicate, flag it
			if ($curStart==$starts[$i] && $curEnd==$ends[$i]) {
				$flag=1;
				last;
			}
			$i++;	#increment counter
		} #end for loop checking for duplicates
		#if not a duplicate, store the information
		if ($flag==0) {
			push(@starts,$curStart);
			push(@ends,$curEnd);
			push(@strands,$curStrand);
		}
	} else { #if current features is not a CDS, store it as such
		push(@otherFeatStarts,$curStart);
		push(@otherFeatEnds,$curEnd);
	}
} #end for loop to compile list of starts/ends
#update log and user on progress
$nCDSs=@starts;
update("done (".$nCDSs." CDS annotations found)\n",$v,$log);



#***************************
if (1==0) {
#***************************


#########################################
# Find predicted pseudogenes
#########################################
#update log and user on progress
update("\tLooking for predicted pseudogenes...",$v,$log);
#declare new variables to store good targets
my @goodLocs;
my @goodGaps;
#run through the current CDSs
my $i=0;
my $convCnt=0;


my $maxHits=0;


foreach my $entry (@starts) {
	#make sure the loop ends properly
	if ($i==$nCDSs) { last; }

print "\n\n----cur seq: ".$seq->seq()."----\n\n";

	my $blastFactory=Bio::Tools::Run::StandAloneBlastPlus->new();

	#grab the sequence of the current CDS
	my $curCDS=$seq->subseq($starts[$i],$ends[$i]);
	my $curCDSobj=Bio::Seq->new(-id=>'curCDS',-seq => $curCDS);
	#now blast it against the genome
	my $curReport=$blastFactory->bl2seq(	-method => 'blastn',
						-query => $seq,
						-subject => $curCDSobj,
						-outfile => 'report.bls');

	die("Died looping through @starts");


	$blastFactory->cleanup;


	my $x=0;

	my $curResult=$curReport;

	my $nHits=$curResult->num_hits;
	if ($nHits>maxHits) { $maxHits=$nHits; }

	#process blast report
	while ($curHit=$curResult->next_hit) {
		
		$x++;

		print "\n\n------\n";
		print "CDS #=$i\nHit #=$x\n";
		print "Max Hits=$maxHits\n";
	}
	$i++; #increment counter
} #end for loop going through each CDS
#if no candidates found, return null
if (@goodLocs==0) {
	#update log and user on progress
	update("none found!\n",$v,$log);
	#return null
	return '';
} else {
	#update log and user on progress
	update("done (".@goodLocs." pseudogene candidates found)\n",$v,$log);
	##update log and user on progress
	#update("\t\t(".$convCnt." convergent CDS pairs found)\n",$v,$log);
}


#***************************
}
#***************************


#########################################
# Now find all places where there are convergently transcribed genes with a gap
# of length between $minGap and $maxGap
#########################################
#update log and user on progress
update("\tLooking for convergent CDS's matching gap parameters...",$v,$log);

#declare new variables to store good targets
#my @goodLocs;
#my @goodGaps;

#run through the current CDSs
my $i=0;
my $convCnt=0;
foreach my $entry (@starts) {
	#make sure the loop ends properly
	if ($i==$nCDSs) { last; }
	#first check if this has a convergent neighbor
	if ($strands[$i]==1 && $strands[$i+1]==-1) {
		$convCnt++;
		#if so, check gap size (difference between convergent ends)
		my $curGap = $starts[$i+1] - $ends[$i];
		if ($curGap > $minGap && $curGap < $maxGap) {
			#if so, add it to the list
			push(@goodLocs,$ends[$i]);
			push(@goodGaps,$curGap);
		} #end gap size check
	} #end convergence check
	$i++; #increment counter
} #end for loop going through each CDS
#if no candidates found, return null
if (@goodLocs==0) {
	#update log and user on progress
	update("none found!\n",$v,$log);
	#return null
	return '';
} else {
	#update log and user on progress
	update("done (".@goodLocs." convergent pairs found)\n",$v,$log);
	#update summary variable for later
	@gapSummaryAllLocs=@goodLocs;
	@gapSummaryAllGaps=@goodGaps;
}

#########################################
# Check if there any other annotations with the gaps
#########################################
#update log and user on progress
update("\tLooking for non-CDS annotations in gap...",$v,$log);
#go through each candidate
my $nCandidates=@goodLocs;
foreach (my $i=$nCandidates-1; $i>=0; $i--) {
	#compute gap boundaries
	my $curGapMin=$goodLocs[$i];
	my $curGapMax=$goodLocs[$i]+$goodGaps[$i];
	#check other features that start in the gap
	my $j=0;
	my $startFlag=0;
	my $endFlag=0;
	foreach (my $j=0; $j<=@otherFeatStarts; $j++) {
		#grab start/end of current feature
		$curStart=$otherFeatStarts[$j];
		$curEnd=$otherFeatEnds[$j];
		#check if the feature starts or ends within the gap
		$startFlag=($curStart > $curGapMin && $curStart < $curGapMax);
		$endFlag=($curEnd > $curGapMin && $curEnd < $curGapMax);
		#if current feature starts or ends within the gap, remove it as a candidate
		if ($startFlag || $endFlag) {
			#remove the current candidate
			splice(@goodLocs,$i,1);
			splice(@goodGaps,$i,1);
			#stop searching for features in this candidate gap
			last;	
		} 
	}
}
#if no candidates remaining, return null
if (@goodGaps==0) {
	#update log and user on progress
	update("no candidates remaining!\n",$v,$log);
	#return null
	return '';
} else {
	#update log and user on progress
	update("done (".@goodGaps." candidates remaining)\n",$v,$log);
	#update summary variable for later
	@gapSummaryOtherAnn=@goodLocs;
#	my @passedOtherAnnotations=@goodLocs;
}

#########################################
# Clear out any targets that have repetitive structure in the gap region
#########################################
#update log and user on progress
update("\tLooking for local repetitive structure...",$v,$log);
#go through each candidate
my $nCandidates=@goodLocs;
foreach (my $i=$nCandidates-1; $i>=0; $i--) {
	#-------------------------
	# Grab the relevant sequence in which to search for local repetitive 
	# structure
	#-------------------------
	#compute search bounds as center of gap +/- the intended homologous 
	# recombination overlap size
	my $gapCenter=$goodLocs[$i]+floor($goodGaps[$i]/2);
	my $locRepLbound=$gapCenter-$homRecombSize;
	my $locRepUbound=$gapCenter+$homRecombSize;

	#if the lower bound is less than 1, set to 1
	if ($locRepLbound<1) { $locRepLbound=1; }
	#grab the sequence
	my $curGapSeq=$seq->subseq($locRepLbound,$locRepUbound);

	#declare some loop variables
	my $hit=0;
	my $locRegLen=length($curGapSeq);
	my $maxWindowSize=floor($locRegLen/$locRep_Reps);
	if ($maxWindowSize>10) { $maxWindowSize=10; }
	#go through each possible window size
	foreach (	my $curWindowSize=$locRep_Nmer; 
			$curWindowSize<=$maxWindowSize; 
			$curWindowSize++ ) {
		#go through each possible window for current window size
		foreach (	my $curWindowInd=0; 
				$curWindowInd<=$locRegLen-$curWindowSize; 
				$curWindowInd++ ) {
			#grab sequence for current window
			my $curWindow = substr $curGapSeq, $curWindowInd, $curWindowSize;
			#check if current window sequence repeats enough
			if ($curGapSeq=~m/($curWindow){$locRep_Reps,}/i) {
				#remove the current candidate
				splice(@goodLocs,$i,1);
				splice(@goodGaps,$i,1);
				#indicate that a hit is found to end looping
				$hit=1;
				last;
			} 
		} #end for loop for possible windows of current size
		#end window size loop if repetitive sequence already found
		if ($hit==1) { last; }
	} #end for loop for possible window sizes
} #end for loop for candidate sites
#if no candidates remaining, return null
if (@goodLocs==0) {
	#update log and user on progress
	update("no candidates remaining!\n",$v,$log);
	#return null
	return '';
} else {
	#update log and user on progress
	update("done (".@goodLocs." candidates remaining)\n",$v,$log);
	#update summary variable for later
	@gapSummaryLocReps=@goodLocs;
}

#########################################
# Clear out any targets that have large repetitive elements nearby
#########################################
#update log and user on progress
update("\tLooking for large-scale repetitive structure...",$v,$log);
#go through each candidate
my $nCandidates=@goodLocs;
foreach (my $i=$nCandidates-1; $i>=0; $i--) {
	#-------------------------
	# Grab the relevant sequence in which to search for local repetitive 
	# structure
	#-------------------------
	#compute search bounds as center of gap +/- the intended homologous 
	# recombination overlap size
	my $gapCenter=$goodLocs[$i]+floor($goodGaps[$i]/2);
	my $lgRepLbound=$gapCenter-$lgRepBound;
	my $lgRepUbound=$gapCenter+$lgRepBound;
	#if the lower bound is less than 1, set to 1
	if ($lgRepLbound<1) { $lgRepLbound=1; }
	#grab the sequence
	my $curRange=$seq->subseq($lgRepLbound,$lgRepUbound);

	#declare some loop variables
	my $hit=0;
	my $lgRegLen=length($curRange);
	#go through each window of the minimum size (note: checking windows of larger
	#	size is redundant in this case)
	$curWindowSize=$lgRepMinSize;
	#go through each possible window for current window size
	foreach (	my $curWindowInd=0; 
			$curWindowInd<=$lgRegLen-$curWindowSize; 
			$curWindowInd++ ) {
		#grab sequence for current window
		my $curWindow = substr $curRange, $curWindowInd, $curWindowSize;
		#grab remaining sequence
		my $curSubRange = substr $curRange, $curWindowInd+1;
		#check if current window sequence repeats
		if ($curSubRange=~m/($curWindow)/i) {
			#remove the current candidate
			splice(@goodLocs,$i,1);
			splice(@goodGaps,$i,1);
			#indicate that a hit is found to end looping
			$hit=1;
			last;
		} 
	} #end for loop for possible windows of current size
} #end for loop for candidate sites
#if no candidates remaining, return null
if (@goodLocs==0) {
	#update log and user on progress
	update("no candidates remaining!\n",$v,$log);
	#return null
	return '';
} else {
	#update log and user on progress
	update("done (".@goodLocs." candidates remaining)\n",$v,$log);
	@gapSummaryLargeReps=@goodLocs;
}

#***************************
if (1==0) {
#***************************

#########################################
# Remove any candidates with potential ORFs nearby
#########################################
#update log and user on progress
update("\tLooking for nearby unannotated/potential ORF's...",$v,$log);
#go through each candidate
my $nCandidates=@goodLocs;
foreach (my $i=$nCandidates-1; $i>=0; $i--) {
	#-------------------------
	# Grab the relevant sequence in which to search for ORFs
	#-------------------------
	#compute the lower search bound as: the 3' end of the gap (which is what
	# is stored in @goodLocs), minus the size of the gap (stored in @goodGaps)
	# and the size of the buffer you want to scan outside the gap
	my $LboundORF=$goodLocs[$i]-$ORFbuffer;
	#compute the upper bound as: the 3' end of the gap plus the buffer size
	my $UboundORF=$goodLocs[$i]+$goodGaps[$i]+$ORFbuffer;
	#if the lower bound is less than 1, set to 1
	if ($locRepLbound<1) { $locRepLbound=1; }
	#actually grab and store the sequence
	my $curGapSeq=$seq->subseq($LboundORF,$UboundORF);
	my $curGapSeqObj=Bio::Seq->new(-seq => $curGapSeq);

	#-------------------------
	# Now run the ORF finder
	#-------------------------
	#initialize ORF Finder factory
	my $EMBOSSfactory=Bio::Factory::EMBOSS->new();
	my $ORFfactory=$EMBOSSfactory->program('getorf');
	#run ORF finder
	my $ORFs=$ORFfactory->run({	-sequence => $curGapSeqObj,
					-outseq => 'orfResults.out',
					-minsize => $minORFsize,
					-maxsize => $maxORFsize, 
					-find => 3 });



print "\n\n----cur seq: ".$seq->seq()."----\n\n";

	#-------------------------
	# If any ORFs found, remove from list
	#-------------------------
	#Check the output file for contents. If it's empty, there were no ORFs
	# found. If there are contents in the file, then ORFs were found. Thus,
	# if the file has contents, remove from the list of candidates
	if (-s 'orfResults.out') {
		#remove the current candidate
		splice(@goodLocs,$i,1);
		splice(@goodGaps,$i,1);
	}
#	#increment index counter
#	$i++;
}
#if no candidates remaining, return null
if (@goodGaps==0) {
	#update log and user on progress
	update("no candidates remaining!\n",$v,$log);
	#return null
	return '';
} else {
	#update log and user on progress
	update("done (".@goodGaps." candidates remaining)\n",$v,$log);
}

#***************************
}
#***************************
############################
#BLAST gaps and look for hits of minGapIdentityWarnLen and minIdentityWarn - warn as necessary
###########################
#update log and user on progress
update("\tChecking minimum gap identities...\n",$v,$log);
#go through each candidate
my $nCandidates=@goodLocs;
my @curIdentWarn;
foreach (my $i=$nCandidates-1; $i>=0; $i--) {

	$curIdentWarn[$i] = 'Y'; 
	#grab the sequence of the current CDS
        my $curSeq=$seq->subseq($goodLocs[$i],$goodLocs[$i] + $goodGaps[$i]);

	my $curSeqobj=Bio::Seq->new(-id=>'curSeq',-seq => $curSeq);
        #now blast it against the genome
	my $blastFactory=Bio::Tools::Run::StandAloneBlastPlus->new();
        my $curReportt=$blastFactory->bl2seq(-method => 'blastn', -query => $curSeqobj, -subject => $seq);


	$blastFactory->rewind_results;
	while (my $curReport=$blastFactory->next_result())
	{
       		#process blast report
	        while ($curHit=$curReport->next_hit) {
			while(my $curHSP=$curHit->next_hsp){
				$len=$curHSP->length('hit');
				$start=$curHSP->start('hit');
				$end=$curHSP->end('hit');
				$name=$curHit->name;
				$num=$curHit->num_hsps;
				next if($goodLocs[$i]==$start);
				
				if($len >= $minGapIdentityWarnLen)
				{	
					my $pctId=$curHSP->percent_identity;
					if($pctID >= minIdentityWarn / 100)
					
					{
						# Failed minimum identity test
						update("\tWARNING. Match with greater than MinIdentity found.\n", $v, $log);
						$curIdentWarn[$i] = 'N';
					}
				}
			}#hsp loop
		}#hit loop
	}#result loop
	$blastFactory->cleanup;
}#main for loop


##############################################
# Find PAM sites
#############################################


update("\tSearching for PAM sites...\n",$v,$log);

my @PAMlocs;
my $nCandidates=@goodLocs;
my $blastFactory=Bio::Tools::Run::StandAloneBlastPlus->new();
my @guideRNAlocs;
my @revguideRNAlocs;
foreach (my $i=$nCandidates-1; $i>=0; $i--) {

	#Grab the candidate
	my $curSeq=$seq->subseq($goodLocs[$i],$goodLocs[$i] + $goodGaps[$i]);
        my $curSeqobj=Bio::Seq->new(-id=>'curSeq',-seq => $curSeq);
	
	# Search curSeq for the PAM site
	while($curSeq =~ m/[ATCGN]GG/gi)
	{
		my $PAMloc = $goodLocs[$i] + $-[0];
		print "PAMloc:$PAMloc\n";
		my $guideRNAcoord = $PAMloc - $guideRNAlen;
		push (@guideRNAlocs, $guideRNAcoord);
	
	}
	# Search for reverse strand PAM sites
	while($curSeq =~ m/CC[ATCGN]/gi)
	{
		my $PAMloc = $goodLocs[$i] + $-[0];
		print "RevPAMloc:$PAMloc\n";
		my $revguideRNAcoord = $PAMloc + 3;
		push (@revguideRNAlocs, $revguideRNAcoord);
	}
	
}
update("\tStoring guide RNA sequences...\n",$v,$log);

my $pass = 0;


#update log and user on progress
update("\tChecking for maximum guide RNA matches...\n",$v,$log);

open(guideRNA,">guideRNAlist.csv") or die ("Failed to open: $!\n");

print guideRNA "Guide RNA start, Guide RNA end,DNA Strand,Guide RNA sequence (5'->3'),BLAST Hits with PAM, Max identity, Max identity w/PAM\n";

#go through each candidate
my $nCandidates=@guideRNAlocs;
print "\tguideRNAlocs:$nCandidates\n";
foreach (my $i=$nCandidates-1; $i>=0; $i--) {

	#grab the sequence of the current guide RNA
	my $guideRNAend = $guideRNAlocs[$i] + $guideRNAlen-1;
        my $curSeq=$seq->subseq($guideRNAlocs[$i],$guideRNAend);
	#print $curSeq."\n";
	my $curSeqobj=Bio::Seq->new(-id=>'curSeq',-seq => $curSeq);
        #now blast it against the genome
	my $blastFactory=Bio::Tools::Run::StandAloneBlastPlus->new();
        my $curReportt=$blastFactory->bl2seq(-method => 'blastn', -query => $curSeqobj, -subject => $seq, -method_args => [ -task=>'blastn-short']);

	$blastFactory->rewind_results;
	while (my $curReport=$blastFactory->next_result())
	{
       		#process blast report
	        while ($curHit=$curReport->next_hit) {
			my $result = "ACCEPT";
			my $pamPresent = "FALSE";
			my $maxident =0;
			my $maxidentwPAM=0;
			while(my $curHSP=$curHit->next_hsp){
				my $ident=$curHSP->num_identical();
				my $len=$curHSP->length('hit');
				$start=$curHSP->start('hit');
				$num=$curHit->num_hsps;
				if ($num<1) # there should be at least 1 hit to itself.
				{
					update("\t WARNING:NO HITS for this guide RNA sequence: $start\n", $v, $log);
				}
				elsif($num == 1)
				{
					$result = "ACCEPT";
					last;
				}  
                                next if($guideRNAlocs[$i]==$start);

				my $seqSearch = $seq->subseq($guideRNAend,$guideRNAend + 3); 
				
				if ($ident > $maxident )
				{
					$maxident = $ident;
				}

				if($seqSearch =~ m/[ATCGN]GG/gi)
				{	
					if ($ident > $maxidentwPAM )
                                	{
                                        	$maxidentwPAM = $ident;
                                	}
			        
					$pamPresent = "TRUE";
					#update("\t guide RNA with PAM found:$start\n", $v, $log);
					if($ident > $maxGuideRNAmatches)
					{		
						#print "$ident\n";
						$result = "REJECT";
						last;
					}
				}
			}#hsp loop
			#update("\t $result this guide RNA sequence: $start.\n", $v, $log);
			if($result eq "ACCEPT")
			{
				print guideRNA "$guideRNAlocs[$i],$guideRNAend,Forward,$curSeq, $pamPresent, $maxident, $maxidentwPAM\n";
				$pass++;
			}

		}#hit loop
	}#result loop
	$blastFactory->cleanup;
}#main for loop


update("\tChecking for maximum guide RNA matches-reverse...\n",$v,$log);
#go through each candidate
my $nCandidates=@revguideRNAlocs;
print "\trevguideRNAlocs:$nCandidates\n";
foreach (my $i=$nCandidates-1; $i>=0; $i--) {

	#grab the sequence of the current guide RNA
	my $guideRNAend = $revguideRNAlocs[$i] + $guideRNAlen-1;
        my $curSeq=$seq->subseq($revguideRNAlocs[$i],$guideRNAend);

	#print $curSeq."\n";
	my $curSeqobj=Bio::Seq->new(-id=>'curSeq',-seq => $curSeq);
        #now blast it against the genome
	my $blastFactory=Bio::Tools::Run::StandAloneBlastPlus->new();
        my $curReportt=$blastFactory->bl2seq(-method => 'blastn', -query => $curSeqobj, -subject => $seq, -method_args => [ -task=>'blastn-short']);

	$blastFactory->rewind_results;
	while (my $curReport=$blastFactory->next_result())
	{
       		#process blast report
	        while ($curHit=$curReport->next_hit) {
			my $result = "ACCEPT";
			my $pamPresent = "FALSE";
			my $maxident =0;
			my $maxidentwPAM=0;
			while(my $curHSP=$curHit->next_hsp){
				my $ident=$curHSP->num_identical();
				my $len=$curHSP->length('hit');
				$start=$curHSP->start('hit');
				$num=$curHit->num_hsps;
				if ($num<1) # there should be at least 1 hit to itself.
				{
					update("\t WARNING:NO HITS for this guide RNA sequence: $start\n", $v, $log);
				}
				elsif($num == 1)
				{
					$result = "ACCEPT";
					last;
				}  
                                next if($revguideRNAlocs[$i]==$start);

				my $seqSearch = $seq->subseq($revguideRNAlocs[$i] - 3,$revguideRNAlocs[$i]); 
				if($ident > $maxident)
				{
					$maxident = $ident;
				}
        

				if($seqSearch =~ m/CC[ATCGN]/gi)
				{	
					if ($ident > $maxidentwPAM )
                                        {
                                                $maxidentwPAM = $ident;
                                        }
					$pamPresent = "TRUE";
					#update("\t guide RNA with PAM found:$start\n", $v, $log);
					if($ident > $maxGuideRNAmatches)
					{		
						#print "$ident\n";
						$result = "REJECT";
						last;
					}
				}
			}#hsp loop
			#update("\t $result this guide RNA sequence: $start.\n", $v, $log);
			if($result eq "ACCEPT")
			{
				my $revcomp = reverse $curSeq;
				$revcomp =~ tr/ATGCatgc/TACGtacg/;
				print guideRNA "$guideRNAend,$revguideRNAlocs[$i],Reverse,$revcomp, $pamPresent, $maxident, $maxidentwPAM\n";
				$pass++;
			}


		}#hit loop
	}#result loop
	$blastFactory->cleanup;
}#main for loop

close(guideRNA);
print "TOTAL PASS:$pass\n";

#############################################
# Wrap things up and write info on good candidates to file
#############################################
#update log and user on progress, as necessary
update(@goodLocs." Candidates passed!\n\n",$v,$log);
#Update targetList.fa file
open(curTargetList,">>targetList.fa") or die ("Failed to open: $!\n");
for (my $i=0; $i<@goodLocs; $i++) {
	#grab sequence +/- gap center for insertion
	my $gapCenter=$goodLocs[$i]-floor($goodGaps[$i]/2);
	my $leftTarg=$gapCenter-$homRecombSize;
	my $rightTarg=$gapCenter+$homRecombSize;
	#if the lower bound is less than 1, set to 1
	if ($leftTarg<1) { $locRepLbound=1; }
	#grab the sequences
	my $leftTargSeq=$seq->subseq($leftTarg,$gapCenter);
	my $rightTargSeq=$seq->subseq($gapCenter+1,$rightTarg);
	print curTargetList "> $genomeFile, location = ".$goodLocs[$i].
		", left arm\n".$leftTargSeq."\n";
	print curTargetList "> $genomeFile, location = ".$goodLocs[$i].
		", right arm\n".$rightTargSeq."\n";
}
close(curTargetList);

#Update targetSummary.txt file
open(targetSummary,">targetSummary.txt") or die ("Failed to open: $!\n");

print targetSummary "Genome \tStart \tStop \tGap Size \tConvergent Gap?".
	"\tPassed Other Annotation Test? \tPassed Local Repetition Test? \tPassed Global Repetition Test? \tPassed Min Gap Identity Test?\n";

for (my $i=0; $i<@gapSummaryAllLocs; $i++) {
	#compute gap info
	my $curStart=$gapSummaryAllLocs[$i];
	my $curStop=$gapSummaryAllLocs[$i]+$gapSummaryAllGaps[$i];
	my $curGapSize=$gapSummaryAllGaps[$i];



	#record which checks it passed/failed
	if ($gapSummaryAllLocs[$i] ~~ @gapSummaryLargeReps) {
		$check1='Y'; $check2='Y'; $check3='Y'; $check4='Y';
	} elsif ($gapSummaryAllLocs[$i] ~~ @gapSummaryLocReps) {
		$check1='Y'; $check2='Y'; $check3='Y'; $check4='N';
	} elsif ($gapSummaryAllLocs[$i] ~~ @gapSummaryOtherAnn) {
		$check1='Y'; $check2='Y'; $check3='N'; $check4='';
	} else {
		$check1='Y'; $check2='N'; $check3=''; $check4='';
	}
	if ($gapSummaryAllLocs[$i] ~~ @goodLocs)
	{
		
		my @indexes= grep { $goodLocs[$_] eq $gapSummaryAllLocs[$i]} 0..$#goodLocs;
		$index=@indexes[0];
		$check5= $curIdentWarn[$index];
	}
	else
	{
		$check5 = '';
	}
	#convert checks array to prinatble format

    my $start=$curStart+1;
	my $stop=$curStop-1;

	print targetSummary "$genomeFile\t$start\t$stop\t$curGapSize\t".
		"$check1\t$check2\t$check3\t$check4\t$check5\n";
	

}
close(targetSummary);

} #end findTargets subroutine
#########################################
# Subroutine update: Updates log file and user on progress, depending on 
#  options chosen by the user via the parameters.txt file.
# Inputs:
#   - update: string containing the update text
#   - v: whether or not to output updates to screen
#   - log: whether or not to output updates to log file
# Outputs:
#   - none
#########################################
sub update {	
	#interpret inputs
	my $update=@_[0];	#string containing update text
	my $v=@_[1];		#output to screen flag
	my $log=@_[2];		#output to log file flag
	#update log and user on progress, as necessary
	if ($log==1) { print logFile $update; }
	if ($v==1) { print $update; }
} #end update subroutine
1;
