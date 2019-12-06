#!/usr/bin/python
from __future__ import print_function #allows colored text?
# This script takes BLAST results and returns back name, frame and amino acid sequence.
import sys
import re
import warnings
import argparse 
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

########################
#   COMMAND OPTIONS    #
########################

parser = argparse.ArgumentParser()
parser.add_argument("in_file", type=file, help="BLAST result file")
parser.add_argument("-o", "--output",  type=file, help="name of output file.")
parser.add_argument("-d", "--debug",  action="store_true", help="name of output file.")
args=parser.parse_args()
debug = args.debug

#Determine outfile name
#input_name = args.in_file.name
#if input_name.find('.'):
#  output_name = input_name.replace


#################
#   FUNCTIONS   #
#################

def store_new_protein(name, frame, expect, score, sequence):
  protein=""">{} | frame={} score={} expect={}
sequence
""".format(name, frame, score, expect, sequence)
  if args.output: 
    with open(args.output, 'w') as out_file:
      out_file.print(protein)
      out_file.close()
  else: print(protein)

#######################
#   REGEXP MATCHING   #
#######################

get_name = re.compile('Query= (.*)\n')
get_score = re.compile(".*Score = (.*) bits.*Expect = (.*),")
get_expect = re.compile('Expect = ([^,]*),')
get_frame = re.compile('.*Frame = ([\+\-][123])')
get_sequence = re.compile('.*Query\s+([0-9]+)\s+([^\s]+)\s+([0-9]+)')
name = frame = expect = score = sequence = ''

state=-1
file_object=args.in_file
for line in file_object: 
  if line.find("Query=")>-1:
    if state==2 or state==1: 
      warnings.warn("unexpected state! ({})".format(state))
    if state==3:
      store_new_protein(name, frame, expect, score, sequence)
    name = frame = expect = score = sequence = ''
    matchObject=re.match(get_name, line)
    name=matchObject.group(1)
    state=1 #"Name gotten"
    file_object.next()

  if state==1 and line.find("Expect =")>-1:  
    matchObject=re.match(get_score, line)
    score=matchObject.group(1)
    expect=matchObject.group(2)
    state=2 #"Expect gotten"
    file_object.next()     

  if state==2 and line.find("Frame")>-1: 
    matchObject=re.match(get_frame, line)
    frame = matchObject.group(1)
    state=3 #"Frame gotten"
    file_object.next()

  if state==3 and line.find("Query ")>-1: 
    if debug: warnings.warn("Line: '{}'".format(line))
    matchObject=re.match(get_sequence, line)
    sequence += matchObject.group(2)

  if debug: 
    print("Current line: '{}', state='{}', name='{}'".format(line, state, name).rstrip())

if state==3:
  store_new_protein(name, frame, expect, score, sequence)

file_object.close()



