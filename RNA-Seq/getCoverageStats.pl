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
# sorted and converted to a pileup file using Samtools. The output is the minimum,
# maximum, mean and median depth of coverage over the entire genome.
#
#!/usr/bin/perl
use strict;
my $pileupfilename=$ARGV[0];

open(PILEUP, $pileupfilename) or die ("Could not open pileup file");
my $min = 999999999999999999999999999;
my $max = -1;
my $mean = 0;
my $median = -1;
my $counter = 0;

use Tie::Array::Sorted::Lazy;

tie my @sortedArray, "Tie::Array::Sorted::Lazy", sub { $_[0] <=> $_[1] };



foreach my $line (<PILEUP>)
{
		my @pileupArray= split '\t', $line;
		my $value= $pileupArray[3];
		if ($value > $max)
		{
			$max = $value;
		}
		if ($value < $min)
		{	
			$min = $value;
		}
		$mean = $mean + $value;
		push @sortedArray, $value;
		$counter ++;
		if ($counter % 200000 == 0)
		{
			print STDERR "Processed 200000 bases\n";
		}
}

$mean = $mean / $counter;
my $length = scalar(@sortedArray);
#print $length;
my $index = 0;
if ($length %2 == 0)
{
	$index = $length /2;
	$median = ($sortedArray[$index] + $sortedArray[$index-1])/2;
}
else
{
	$index = int($length /2);
	$median = $sortedArray[$index];
}

print "Min: $min, Max: $max, Mean: $mean, Median: $median\n"; 
close(PILEUP);
