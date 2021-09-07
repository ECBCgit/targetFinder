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
#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use POSIX;
use Bio::Seq;
use Bio::SeqIO;
#########################################
# Pull in other perl scripts
#########################################
require 'Code/find_targets.pl';
require 'Code/extract_params.pl';
#########################################
# Load input information
#########################################
# Change directory to correct project folder
my $projectFolder=$ARGV[0];
chdir("Projects/".$projectFolder);
my $projPath=getcwd();
#----------------------------------------
# Load parameters from parameter file
#----------------------------------------
my @params=extractParams();
# Declare parameter groups
my @initParams=();
my @genParams=();
my @targetParams=();
# Fill in each parameter group
my $groupCnt=0; my $paramCnt=0;
foreach my $param (@params) {
	if ($param=~/x/) { $groupCnt++; $paramCnt=0 }
	elsif ($groupCnt==0) { $initParams[$paramCnt]=$param; $paramCnt++; }
	elsif ($groupCnt==1) { $genParams[$paramCnt]=$param; $paramCnt++; }
	elsif ($groupCnt==2) { $targetParams[$paramCnt]=$param; $paramCnt++; }
}
#----------------------------------------
# Set up logs
#----------------------------------------
chdir($projPath."/Logs") or die "Failed to open directory Logs: $!";
# If instructed to, remove old log files
my $clearLogsFlag=$initParams[0];
if ($clearLogsFlag==1) { 
	foreach my $logFile (<*>) 
		{ unlink($logFile); 
	} 
}
# Grab system time
(my $sec,my $min,my $hour,my $mday,my $mon,my $year,my $wday,my $yday,
	my $isdst)=localtime(time);
# Format timestamp
my $timeStamp=sprintf "%4d_%02d_%02d_%02d_%02d_%02d",
		$year+1900,$mon+1,$mday,$hour,$min,$sec;
# Open log file
open(logFile,">".$timeStamp.".log") or die "Failed to open results file: $!";
chdir($projPath);

#########################################
#########################################
my $dir=getcwd();
chdir('targetGenomes');
my @genomes=<*>;
chdir($dir);

#########################################
# Find target locations for each of the target genomes.
#########################################
#Go through each target genome and find potential targets
my @targets;
foreach my $genome (@genomes) {
	#skip files ending in ~
	if ($genome=~/\w$/) {
		# Find and store targets
		findTargets($genome,@genParams,@targetParams);
	}
}
#########################################
# Wrap everything up
#########################################
print logFile "\nScript end.\n\n";
print "\nScript end.\n\n";
close(logFile);


