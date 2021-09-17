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
# Introduction
#########################################
# This document describes how to use targetFinder. The code will
# generate a list of putatively neutral insertion sites in prokaryotes, 
# along with a list of guide RNA sequences for the gap. 
#
# Running the code
#
# This software was tested on CentOS 7, and was tested on machines with at
# least 32 GB of RAM. However, as no large BLAST databases are used, it 
# should run with significantly lower amounts of RAM. A modern CPU is 
# recommended, however only the BLAST searches utilize multi-threaded
# processing. The default behavior is to use all cores available for BLAST. 
#
# Dependencies:
#
# * perl 5.8 or higher
# * Bioperl. Information on installing Bioperl can be found at
#   https://bioperl.org/INSTALL.html Specific Bioperl modules used in this program
#   include:
#	- StandAloneBlastPlus (bl2seq, blastn)
#	- EMBOSS
#
# Test data is located in the Projects directory.	
#
#########################################
# User Guide
#########################################
# The first step is to create a new project directory. This directory needs to 
# 	be placed under the Barcode/Projects/ directory. There is a directorie
#	there called EcoliMG1655 as an example.
#
# This folder should contain the following:
#	(1) Subdirectory called targetGenomes containing the target genome. 
#		This must be in genbank format.
#	(2) Text file containing the script parameters, called "parameters.txt".
#		This file should be copied from the EcoliMG1655  
#		project and tweaked as appropriate. Each parameter is described 
#		in the "parameters.txt" file in the project directory.
#	(3) Subdirectory called Logs, which may start as empty.
#
# Usage: 
# 	from the targetFinder/ directory, type:
#
#	perl targetFinder.pl "project_name" 
#
#	Where:
#	- project_name = name of the subdirectory of the current project, located
#		at Barcode/Projects/project_name
#	
# Outputs:
#	- targetSummary.txt: list of candidate gaps found, along with whether each
#		of the tests passed for each gap.
#	- targetList.fa: fasta file that contains the left and right arms of gaps that
#		passed all of the tests.
#	- guideRNAlist.csv: csv file containing guide RNA sequences that were found in 
#		each of the gaps.
#	- log files: complete description of what the script did (see
#	 	parameters.txt for some options)
#
#########################################
