#! /usr/bin/env python2.7

#Code produced by Chris Mutzel
#Christopher.Mutzel@Tufts.edu
#Tufts University, Medford, MA 02155
#Last updated January 27, 2011

#Modules
import os
import csv 

work_dir = "/Users/christophermutzel/Documents/School/Power Systems/";
in_file = "30BusInputData.csv";

#Make output file name
out_file = work_dir + in_file + "_matrix";

#Open the file and read in as csv

with open(work_dir + in_file, 'rb') as f:
    reader = csv.reader(f)
    for row in reader:
        print row
        
#End of line        