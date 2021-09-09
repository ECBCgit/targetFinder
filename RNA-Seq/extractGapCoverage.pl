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
###
#
# This perl script takes a pileup file that has been created by Bowtie2, and then
# sorted and converted to a pileup file using Samtools. It also requires a gap file
# from targetFinder. The output is a file with the gap range, and the depth of coverage
# at each base position, in a comma separated format.
#
#!/usr/bin/perl
use strict;
my $gapfilename=$ARGV[0];
my $pileupfilename=$ARGV[1];

open(GAPS, $gapfilename) or die ("Could not open gaps file");
open(PILEUP, $pileupfilename) or die ("Could not open pileup file");

print STDERR "Loading pileup...\n";
	
my @pileup=<PILEUP>;
my $flank=1000;


foreach my $line (<GAPS>)
{
	my @array =split '\t', $line;
	print "Gap: $array[1] $array[2]\n";

	my $start=-1;
	if ($array[1]-$flank >=0)
	{
		$start = $array[1]-$flank;
	}
	else
	{
		$start = 0;
	}
	my $end=-1;
	if ($array[2]+$flank < scalar(@pileup))
	{
		$end = $array[2]+$flank;
	}
	else
	{
		$end = scalar(@pileup) -1;
	}

	
	for (my $i=$start;$i <= $end; $i++)
        {
		print "$i";
		if ($i < $end)
		{	
			print ",";
		}
	}
	print "\n";
	for (my $i=$start; $i <= $end; $i++)
	{
		my $pileupLine=$pileup[$i-1];
		my @pileupArray= split '\t', $pileupLine;
		print "$pileupArray[3]";
		if ($i < $end)
		{	
			print ",";
		}
	}
	print "\n\n";
}
close(PILEUP);
close (GAPS);
