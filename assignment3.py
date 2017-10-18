#!/usr/bin/env python
#import numpy & be able to call it as "np"
import numpy as np
import scipy as sp
import scipy.stats

#import sys
import sys

# turn on debug mode
Debug = True

# helper function that calculates confidence interval and ranges
def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h

#import autocorrelation.py script provided; includes numpy; allows acf command
#import autocorrelation.py
import autocorrelation

#explanation of what this code does if you make an input error by
#having less than three arguments
if len(sys.argv) < 3:
   print("")
   print("My test script to calculate autocorrelation scores for amino acid fasta files.")
   print("Usage: args.py -i <inputfile> -o <outputfile> -w <window size>")
   print("-i: input file")
   print("-o: output file")
   print("")
   sys.exit()

# creating infile and outfile, enable usage of arguments
for i in range(len(sys.argv)):
    # input file arg
    if sys.argv[i] == "-i":
        # store input file name
        InFileName = sys.argv[i+1]
    # output file arg
    elif sys.argv[i] == "-o":
        # store output fle name
        OutFileName = sys.argv[i+1]
    # window size arg
    elif sys.argv[i] == "-w":
        # store window size
        window = sys.argv[i+1]


# set file name
HydropathyFileName = "amino_acid_hydropathy_values.txt"
# open hydropathy values file
HydropathyFile = open(HydropathyFileName, 'r')
# set data array
Data=[]
# set Hydropathy dictionary
Hydropathy={}

# initialize line number at 0
LineNumber = 0

# loop thru Hydropathy file and associate value with amino acid
for Line in HydropathyFile:
    if(LineNumber>0):
        Line = Line.strip("\n")
        Data = Line.split(",")
        Hydropathy[Data[1]]=float(Data[2])
    LineNumber = LineNumber + 1
HydropathyFile.close()

#open the infile and outfile for our data
InFile = open(InFileName, 'r')

# initialize line number (1 to make it more human-readable)
lineNumber = 1
# initialize ProtName
ProtName = ""

# loop thru the lines to get the sequence
for line in InFile:
  if(lineNumber == 1):
    ProtName = line.strip('>')
    ProtName = ProtName.strip()

  # get lines after first line
  if(lineNumber > 1): 
    ProtSeq = line.strip('\n')
  
  # increment
  lineNumber = lineNumber + 1

# close the included file
InFile.close()

# open output file in Append mode
OutFile = open(OutFileName,'a')


# creating window size as raw input, window size we will put in when we run the script in UNIX is 9
window=int(window)
Value=0
window_counter=0

HydropathyValues = []

# window
for i in range(len(ProtSeq)):
    Value+=Hydropathy[ProtSeq[i]]
    if(i>(window-1) and i<=(len(ProtSeq)-window)):
        # value for the window
        Value=Value-Hydropathy[ProtSeq[i-window]]
        #OutString = "%d,%.2f" % (window_counter, Value)
        #OutFile.write(OutString + "\n")

    HydropathyValues.append(Value)
    window_counter+=1

print(HydropathyValues)

# get ACF values from hydropathy values
acf = list(autocorrelation.acf(HydropathyValues))

# calculate the length of acf list for confidence interval (95%)
acfLength = len(acf)

# list of True Mean, and range of 95% confidence interval values
acfMeanAndRanges = mean_confidence_interval(acf)

# set count at 0
count = 0

# loop through each value in acf list
for i in acf:
  # check if i is between the values
  if acfMeanAndRanges[1] <= i <= acfMeanAndRanges[2]:
    True #arbitrary true
  else:
    # add to count b/c outside confidence interval
    count += 1

# calculate proportion of acf values outside confidence interval
proportion = count/acfLength

# tell user how many outside confidence interval for this particular file
print(str(proportion) + " of correlation values are outside of the 95% confidence interval.")

# initialize var
outOfConfidenceInterval = "No"

# check if proportion of values is outside of confidence interval
if proportion > 0.05:
  outOfConfidenceInterval = "Yes"

# output protein name, proportion acf values outside confidence interval, yes/no, and max hydropathy val
OutString = "%s,%.2f, %s, %.2f" % (ProtName, proportion, outOfConfidenceInterval, max(HydropathyValues))
# write to output file, add line break
OutFile.write(OutString + "\n")

# close output file
OutFile.close()
