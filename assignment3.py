#!/usr/bin/env python
#import numpy & be able to call it as "np"
import numpy as np

#import sys
import sys

# turn on debug mode
Debug = True

#import autocorrelation.py script provided; includes numpy; allows acf command
#import autocorrelation.py
import autocorrelation

#explanation of what this code does if you make an input error by
#having less than three arguments
if len(sys.argv) < 2:
   print("")
   print("My test script to calculate autocorrelation scores for amino acid fasta files.")
   print("Usage: args.py -i <inputfile> -o <outputfile>")
   print("-i: input file")
   print("-o: output file")
   print("")
   sys.exit()

#creating infile and outfile, enable usage of arguments
for i in range(len(sys.argv)):
    if sys.argv[i] == "-i":
        InFileName = sys.argv[i+1]
    elif sys.argv[i] == "-o":
        OutFileName = sys.argv[i+1]


# open hydropathy values
HydropathyFileName = "amino_acid_hydropathy_values.txt"
HydropathyFile = open(HydropathyFileName, 'r')
Data=[]
Hydropathy={}
LineNumber = 0

for Line in HydropathyFile:
    if(LineNumber>0):
        Line = Line.strip("\n")
        Data = Line.split(",")
        Hydropathy[Data[1]]=float(Data[2])
    LineNumber = LineNumber + 1
HydropathyFile.close()

print(Hydropathy)


#open the infile and outfile for our data
InFile = open(InFileName, 'r')

lineNumber = 1

# loop thru the lines to get the sequence
for line in InFile:

  # get lines after first line
  if(lineNumber > 1): 
    ProtSeq = line.strip('\n')
  
  lineNumber = lineNumber + 1

InFile.close()

print(ProtSeq)
#print(autocorrelation.acf(list(sequence).astype(np.float)))

OutFile = open(OutFileName,'w')


# creating window size as raw input, window size we will put in when we run the script in UNIX is 9
window = input("Window size?")
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
        OutString = "%d,%.2f" % (window_counter, Value)
        OutFile.write(OutString + "\n")

    HydropathyValues.append(Value)
    window_counter+=1


acf = list(autocorrelation.acf(HydropathyValues))

acfLength = len(acf)

acfSquare = np.sqrt(acfLength)

count = 0

print(acf)
# loop through each 
for i in acf:
  if i > acfSquare or i < acfSquare:
    count = count + 1

print(count)


InFile.close()
OutFile.close()
