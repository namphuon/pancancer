from __future__ import division
import pancancer
import sys, os, subprocess
import hg19util as hg19
import json,requests
import argparse
import tempfile,shutil
import glob
import pysam
import os,random,tempfile,shutil
import os.path
import copy,math,bisect
from sets import Set
from collections import Counter
from argparse import ArgumentParser
from googleapiclient.discovery import build
import httplib2
import pprint

chrs = ["%d" % x for x in xrange(1,23)]
chrs.extend(['Y','X'])

def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument('-b', '--bam_file', default="",
                      help='''bam file''')
  parser.add_argument('-s', '--segment_file', default="",
                      help='''segment file''')
  parser.add_argument('-o', '--output_file', default="",
                      help='''output file''')
  parser.add_argument('-f', '--flagstat_file', default="",
                      help='''flagstat file''')                      
  options = parser.parse_args()
  return options
  
def compute_depth(bam_file, intervals, output_file, step_size = 100):
  bam = pysam.Samfile(bam_file, 'rb')
  chr_prefix = bam.references[0].find('chr') != -1
  output = open(output_file, 'w')  
  for i in intervals:    
    curr_pos_idx = 0
    curr_depth = 0
    prev_pos = -1
    for pileupcolumn in bam.pileup(i.chrom if chr_prefix else i.chrom.replace('chr',''), i.start, i.end,stepper='all', truncate=True):
      if prev_pos == -1:
        prev_pos = pileupcolumn.pos
      curr_depth+=pileupcolumn.n
      if pileupcolumn.pos - prev_pos >= step_size:
        output.write("%d,%f\n" % (int((pileupcolumn.pos+prev_pos)/2), curr_depth / (pileupcolumn.pos - prev_pos)))
        prev_pos = pileupcolumn.pos
        curr_depth = 0
  output.close()

def parse_bed(segment_file):
  input = open(segment_file, 'r')
  amplicons = hg19.interval_list()
  for line in input:
    res = line.strip().split('\t')
    amplicons.append(hg19.interval(res[0], int(res[1]), int(res[2])))
  return amplicons


def read_count(bam_file, amplicons, output_file):
  bam = pysam.Samfile(bam_file, 'rb')
  prefix = pancancer.prefix_chr(bam)
  output = open(output_file,'w')
  for a in amplicons:
    counted = {}
    idx = 0
    for read in bam.fetch(a.chrom if prefix else a.chrom.replace('chr',""), a.start, a.end):
      if not read.is_unmapped:
        foo = counted.setdefault(read.qname,"")
        idx+=1
    output.write('%s,%s:%d-%d,%s,%d,%d\n' %
                (a.info['sample'],a.chrom if prefix else a.chrom.replace('chr',""),a.start,a.end,a.info['id'],len(counted.keys()),idx))    
  output.close()

if __name__ == '__main__': 
  options = parse_args()
  #amplicons = pancancer.parse_segment_file(options.segment_file)
  amplicons = parse_bed(options.segment_file)
  #if options.flagstat_file != '':
  #  pancancer.samtools_flagstat(options.bam_file, options.flagstat_file, threads=1)  
  #read_count(options.bam_file, amplicons, options.output_file)
  compute_depth(options.bam_file, amplicons, options.output_file)
  #sys.exit(0)