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
#!/usr/local/bin/perl
use strict;
use warnings;

print "";

sub extractParams {

my @initParams=();
my @genParams=();
my @targetParams=();

open(paramFile,"parameters.txt") or die "Failed to open parameters file: $!";
while (<paramFile>) {
	my $curLine=$_;
	#First grab initialization parameters
	if ($curLine =~ m/initialize_(\S*) = (\d*)/) {
		my $hit=$2;
		if ($1=~/clearLogs/) { $initParams[0]=$hit; }
		elsif ($1=~/EMBOSSpath/) { system("export PATH=\$PATH:$hit"); }
	#Now grab general parameters
	} elsif ($curLine =~ m/general_(\w*) = (\d*)/) {
		my $hit=$2;
		if ($1=~/verbose/) { $genParams[0]=$hit; }
		elsif ($1=~/logs/) { $genParams[1]=$hit; }	
	#Now grab target parameters
	} elsif ($curLine=~m/target_(\w*) = (\d*)/) {
		my $hit=$2;
		if 	($1=~/minGap/) { 	$targetParams[0]=$hit; }
		elsif 	($1=~/maxGap/) { 	$targetParams[1]=$hit; }
		elsif 	($1=~/homRecombSize/) {	$targetParams[2]=$hit; }
		elsif 	($1=~/locRep_Nmer/) { 	$targetParams[3]=$hit; }
		elsif 	($1=~/locRep_Reps/) { 	$targetParams[4]=$hit; }	
		elsif 	($1=~/ORFbuffer/) { 	$targetParams[5]=$hit; }
		elsif 	($1=~/minORFsize/) { 	$targetParams[6]=$hit; }
		elsif 	($1=~/maxORFsize/) { 	$targetParams[7]=$hit; }
		elsif 	($1=~/lgRepBound/) { 	$targetParams[8]=$hit; }
		elsif 	($1=~/lgRepMinSize/) { 	$targetParams[9]=$hit; }
		elsif 	($1=~/lgRepMaxSize/) { 	$targetParams[10]=$hit; }
		elsif   ($1=~/min_gap_identity_warning_length/) {  $targetParams[11]=$hit; }                
		elsif   ($1=~/min_identity_warning/) {  	   $targetParams[12]=$hit; }                
		elsif   ($1=~/pam_site_sequence/) {    		   $curLine=~m/target_(\w*) = ([ATCGNatcgn]*)/;
	      					                   $targetParams[13]=$2; }
		elsif   ($1=~/guide_RNA_length/) {  		   $targetParams[14]=$hit; }                
		elsif   ($1=~/max_guide_RNA_matches/) {  	   $targetParams[15]=$hit; }                
		elsif   ($1=~/small_RNA_seq_window/) {  	   $targetParams[16]=$hit; }                
		elsif   ($1=~/large_RNA_seq_window/) {  	   $targetParams[17]=$hit; }                
	}
}

#Compile into single output w/ junk character spacers for later parsing
my @output=(@initParams,'x',@genParams,'x',@targetParams);
}
