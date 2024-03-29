#!/usr/bin/python
from __future__ import print_function #allows colored text

# This script runs bowtie2 with different settings, doing a random-walk towards the best settings for assembling a chromosome from a fastq file.
# Files are from sequence capture in malaria.
# Bowtie2 matches the reads to chromosome Ht0017, according to slightly mutated settings.
# The results are compared to the "correct" sequence, - reads_correct.fastq - produced in Geneious.
# If the result from the new setting is better than the previous run, the new settings are saved as the template.
# Otherwise, the script loads the old settings and try another random mutation.
# The script halts after a predetermined number of steps.
# All settings and results, including those from previous runs, are saved in settings.txt 

import sys
import re
import os
import warnings
import argparse 
import time
from random import randint
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

#Skip traceback and just print the warning
def custom_formatwarning(msg, *a): 
  if debug:
    print("\nwarning: " + str(msg) + '\n')
  return "warning: " + str(msg) + '\n'
warnings.formatwarning = custom_formatwarning

'''
Typical settings: 
-D 20 -R 3 -N 1 -L 20 -i S,1,0.50
so: 
20 3 1 20 0.50 [match] [score] [time] 
'''

##################
### VARIABLES ####
##################

resolution=0.5
pipeline="/home/roland/bin/pipeline"
settings_file="settings.txt"
settings_dictionary={}
target_fastq="merged_Ht0017_reads.fastq"
target_ref="HtScaffold0001:67143-68859.fasta"
target_index="HtIndex"
target_output="scaffolds/consensus"
target_comparison="Ht0017.fasta"
target_sample="512022_S1_L001"
basic_settings="--local -N 1 -q --threads 3"
nr_of_tries=1000
previous_setting="20 2.2 1 23.2 0.28 0.04"
previous_score=0
current_setting=""
current_score=0
setting=""
folder="."
verbose=0
D_step=1*resolution
R_step=0.5*resolution
L_step=1*resolution
I1_step=0.1*resolution
I2_step=0.1*resolution

######################
### POPULATE ARRAY ###
######################

#Open settings file (or create it)
if (os.path.isfile(settings_file)):
  file_object=open(settings_file, "r")
  for line in file_object: #For every line, create entries in array 
    (setting_D, setting_R, setting_N, setting_L, setting_I1, setting_I2, score)=line.split()
    name="{} {} {} {} {} {}".format(setting_D, setting_R, setting_N, setting_L, setting_I1, setting_I2)
    settings_dictionary[name]=score
    if score>previous_score: 
      previous_setting=name
      previous_score=score
  file_object.close()
file_object=open(settings_file, "a")
print("best setting ({}) has score {}".format(previous_setting, previous_score))


#################
### FUNCTIONS ###
#################


################################
### function to mutate a value
def mutate(value, minimum, maximum, step): 
  dice=randint(1, 3)
  change=step*randint(1, 3)
  if (dice==1): #1/3 chance, increase value 
    if (value+change)<=maximum: 
      return value + change
    return maximum
  elif (dice==2): #1/3 chance, decrease value
    if (value-change)>=minimum: 
      return value - change
    return minimum
  return value # 1/3 chance, do nothing   


################################
### function to evaluate result 
def calculate_score(consensus_file, execution_time):
  global target_comparison
  out_file="scaffolds/temp"

  #run stretcher (compare to 'correct' sequence) 
  command="[[ -e {0} ]] && rm {0}; stretcher {1} {2} {0}".format(out_file, consensus_file, target_comparison)
  os.system( command )

  #get [number of bases matching], and [mismatches], put into variables
  command="cat scaffolds/temp | grep \"# Score\" | grep \"[1-9]*\" -o > score.txt"
  os.system( command )
  file_object=open("score.txt", "r")
  score=int(file_object.read()) - 10*execution_time
  if (verbose): print("score = '{}'".format(score))

  # TO BE DONE: 
  # set penalty according to mismatches and execution time
  # (length, mismatches)=() 
  #score=(length - 3*mismatches - execution_time/300)

  return score

#############################
### function to run bowtie
def run_bowtie2(settings, current_settings):
  if (verbose): print("run_bowtie2, target_output={}".format(target_output))
  #if these settings have not already been tried, return score directly
  if (settings in settings_dictionary): 
    return settings_dictionary[settings]

  
  #run bowtie2 
  if (verbose):print("Running bowtie with setting '{}'".format(current_settings))
  start_time=time.time()
  pre_command="bowtie2 {} {} -U {} -x {} -S {}.sam ".format(basic_settings, settings, target_fastq, target_index, target_output)
  if (verbose==0): command= pre_command+" > /dev/null 2> /dev/null"
  if (verbose):print("command="+command)
  os.system( command )
  execution_time=(start_time-time.time())

  #sort and index SAM 
  if (verbose): print("indexing SAM-file")
  command="bash {3}/SAM_to_indexed_BAM.sh {0}.sam {0}.bam {1} {2}".format(target_output, target_fastq, target_sample, pipeline)
  if (verbose): print("command={}".format(command))
  os.system( command )

  #make consensus
  if (verbose): print("making consensus")
  command="bash {}/consensus_from_samfile.sh {}.bam {} {}.fasta".format(pipeline, target_output, target_comparison, target_output)
  if (verbose): print("command="+command)  
  os.system( command )

  #evaluate result
  if (verbose): print("evaluating result")
  score=calculate_score("{}.fasta".format(target_output), execution_time)
  
  #store setting and its result. Then return score
  file_object.write("{} {}\n".format(current_settings, score))
  return score
  
###########################
### function to run match 
def run_matcher(iteration):
  global current_setting, previous_setting, current_score, previous_score, setting
  #generate settings based on modifying previous settings.

  #Get previous setting
  #print("Previous setting is '{}'".format(previous_setting))
  (setting_D, setting_R, setting_N, setting_L, setting_I1, setting_I2)=map(float, previous_setting.split())
  
  #Every tenth try, take a bigger leap, to get out of local maximums.
  step=1  
  if (iteration%10==9): 
    step=5
    print("\nTaking bigger step...")
  elif (iteration%5==4): 
    step=0.2
    print("\nTaking smaller step...")  

  setting_D = mutate(setting_D, 0, 30, D_step*step)
  setting_R = mutate(setting_R, 0, 5, R_step*step)
  setting_L = mutate(setting_L, 0, 30, L_step*step)
  setting_I1 = mutate(setting_I1, 0, 5, I1_step*step)
  setting_I2 = mutate(setting_I2, 0, 3, I2_step*step)
  if randint(0,5)==5:
    if setting_N==1: 
      setting_N=0
    else: setting_N=1

  #Every twentienth try, average setting values with one of the presets, to keep on track.
  #TBD

  #Delete output folder 
  command="rm -r scaffolds; mkdir scaffolds"
  os.system( command )
  
  #Run bowtie2 and get result
  setting = "-D {} -R {} -N {} -L {} -i S,{},{}".format(setting_D, setting_R, setting_N, setting_L, setting_I1, setting_I2)
  current_setting= "{} {} {} {} {} {}".format(setting_D, setting_R, setting_N, setting_L, setting_I1, setting_I2)
  current_score=run_bowtie2(setting, current_setting)
  
  #if current result is worse than previous result, switch back.
  if (current_score>previous_score):
    previous_setting=current_setting
    previous_score=current_score
    print("\n*** Found new best setting: {} has score {}. ***\n".format(setting, current_score))
  else:  
    print("Rejected setting: {} has score {}, less than {}.".format(setting, current_score, previous_score))




#################
### MAIN LOOP ###
#################

#For up to nr_of_times tries, run matcher
for x in range(1, nr_of_tries):
  run_matcher(x)
print("\nbest setting = {} at score {}".format(previous_setting, previous_score))
