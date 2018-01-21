from __future__ import division
import pickle
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


PANCANCER_DIR=os.environ['PANCANCER']
TEMP_DIR= os.environ['TEMP_DIR'] if 'TEMP_DIR' in os.environ else None
bam_json_file = "%s/metadata/hnsc.bam.json" % PANCANCER_DIR
cnv_json_file = "%s/metadata/hnsc.masked_cnv.json" % PANCANCER_DIR

# the CLIENT_ID for the ISB-CGC site
CLIENT_ID = '907668440978-0ol0griu70qkeb6k3gnn2vipfa5mgl60.apps.googleusercontent.com'
# The google-specified 'installed application' OAuth pattern
CLIENT_SECRET = 'To_WJH7-1V-TofhNGcEqmEYi'
# The google defined scope for authorization
EMAIL_SCOPE = 'https://www.googleapis.com/auth/userinfo.email'
# where a default credentials file will be stored for use by the endpoints
DEFAULT_STORAGE_FILE = os.path.join(os.path.expanduser("~"), '.isb_credentials')
diseases = [d for d in os.listdir('%s/analyses/tumor' % PANCANCER_DIR) if d.find('TCGA') != -1]
diseases.sort()

def build_bed_files(file = '%s/analyses/breakpoints.csv' % PANCANCER_DIR):
  input = open(file,'r')
  beds = {}
  for line in input:
    res = line.strip().split(',')    
    i1 = hg19.interval(res[1],int(res[2])-10000,int(res[3])+10000)
    i2 = hg19.interval(res[4],int(res[5])-10000,int(res[6])+10000)
    foo = beds.setdefault(res[7],hg19.interval_list())
    beds[res[7]].append(i1)
    beds[res[7]].append(i2)
  for sample in beds.keys():
    hits = len(beds[sample])
    new_hits = -1
    while new_hits != hits:
      beds[sample].sort()
      my_list = hg19.interval_list()    
      founds = []
      for i in beds[sample]:
        if i in founds:
          continue      
        matches = [h[0] for h in beds[sample].intersection([i])]
        if len(matches) == 1:
          my_list.append(matches[0])
        else:
          new_int = hg19.interval(matches[0].chrom, min([m.start for m in matches]), max([m.end for m in matches]))
          founds.extend(matches)
          my_list.append(new_int)
      beds[sample] = my_list
      hits =  len(beds[sample])
      new_hits = len(my_list)
    output = open('%s/analyses/tmp/%s.bed' % (PANCANCER_DIR, sample),'w')
    for i in beds[sample]:
      output.write('%s\t%d\t%d\n' % (i.chrom, i.start, i.end))
    output.close()    

def load_bed_file(bed_file, header = False, map = None):
  input = open(bed_file, 'r')
  interval_list = hg19.interval_list()
  if header == True:
    input.next()
  lines = [line.strip().split('\t') for line in input]
  if map == None:
    foo = [interval_list.append(hg19.interval(line[0],int(line[1]),int(line[2]),info={'data':line})) for line in lines]
  else:
    foo = [interval_list.append(hg19.interval(line[map['chr']],int(line[map['start']]),int(line[map['end']]),info=dict([(k,line[v]) for (k,v) in map.items()]))) for line in lines]    
  interval_list.sort()
  return interval_list

def load_everything():
  ensembl = load_ensemble_genes()
  diseases = {}
  [diseases.setdefault(data_map[k][key_map['WGS__cases__project__project_id']],[]).append(k) for k in data.keys() if k in keeper and keeper[k][key_map['rnaseq_gsc']] != '']
  output = open('/pedigree2/projects/namphuon/programs/pancancer/analyses/intervals.csv','w')
  output.write("Disease,Patient,Amplicon,Interval,IsAmplicon,PartAmplicon,IsCycle,FPKM,CNV\n")
  for disease in diseases:
    #Collect all segments by patient, and their total copy count, and whether or not that segment
    #is on a cycle for an amplicon
    interval_map = {}
    intervals = hg19.interval_list()
    for patient in diseases[disease]:
      for amplicon in data[patient][1].keys():
        if amplicon == 'amplicons':
          continue
        #Now classify the amplicon region
        p_interval = hg19.interval_list()        
        for s in data[patient][1][amplicon]['interval_map']:
          segment = copy.copy(data[patient][1][amplicon]['interval_map'][s])
          segment.info['copy'] = 0
          segment.info['cycle'] = False
          segment.info['patient'] = patient
          segment.info['amplicon'] = amplicon
          segment.info['interval'] = s
          segment.info['segments'] = Set()
          p_interval.append(segment)
        p_interval.sort()
        for cycle in data[patient][1][amplicon]['cycle_map'].keys():
          is_cyclic = data[patient][1][amplicon]['cycle_map'][cycle]['cycle'][0][0] != '0'
          for seg in data[patient][1][amplicon]['cycle_map'][cycle]['cycle']:
            if seg[0:-1] != '0':
              #Find out which intervals it intersects, hopefully only one
              hits = p_interval.intersection([data[patient][1][amplicon]['segment_map'][seg[0:-1]]])
              if len(hits) > 1:
                print "Fucked"
                break    
              elif len(hits) == 0:
                print "Double fucked"
                break
              hits[0][0].info['copy']+=(data[patient][1][amplicon]['cycle_map'][cycle]['copy_count'] *(hits[0][1].end-hits[0][1].start))/ (hits[0][0].end-hits[0][0].start)
              hits[0][0].info['cycle'] = hits[0][0].info['cycle'] or is_cyclic
        intervals.extend(p_interval)
    intervals.sort()          
    for i in intervals:
      hits = [h[0] for h in ensembl.intersection([i])]
      for patient in diseases[disease]:
        #Expression of all transcripts that intersect
        fpkm = sum([data[patient][-1][g.info['gene_id']] for g in hits if g.info['gene_id'] in data[patient][-1]])
        copy_count = [(h[0].end-h[0].start,2**float(h[0].info['seg_mean'])) for h in data[patient][0].intersection([hg19.interval(i.chrom.replace('chr',''), i.start, i.end)])]
        if len(copy_count) > 0:
          copy_count = sum([v[0]*v[1] for v in copy_count])/sum([v[0] for v in copy_count])
        else:
          copy_count = 2
        #Copy count, first check if any intersect
        matches = [h[0] for h in intervals.intersection([i]) if h[0] != i and h[0].info['patient'] == patient]
        #if len(matches) != 0:
        #  copy_count = sum([(v.end-v.start)*v.info['copy'] for v in matches])/sum([(v.end-v.start) for v in matches])
        output.write('%s,%s,%s,%s:%d-%d,%s,%s,%s,%f,%f\n'% (disease,patient,i.info['amplicon'],
                      i.chrom,i.start,i.end,
                     patient==i.info['patient'], patient==i.info['patient'] or len(matches) > 0, i.info['cycle'] if i.info['patient']==patient 
                    else len([m for m in matches if m.info['cycle'] == True]) > 0, fpkm, copy_count))
  output.close()

  output = open('/pedigree2/projects/namphuon/programs/pancancer/analyses/segments.csv','w')
  output.write("Disease,Patient,Amplicon,Segment,IsAmplicon,PartAmplicon,IsCycle,FPKM,CNV\n")
  for disease in diseases:
    #Collect all segments by patient, and their total copy count, and whether or not that segment
    #is on a cycle for an amplicon
    interval_map = {}
    intervals = hg19.interval_list()
    for patient in diseases[disease]:
      for amplicon in data[patient][1].keys():
        if amplicon == 'amplicons':
          continue
        #Now classify the amplicon region
        p_interval = hg19.interval_list()        
        for s in data[patient][1][amplicon]['interval_map']:
          segment = copy.copy(data[patient][1][amplicon]['interval_map'][s])
          segment.info['copy'] = 0
          segment.info['cycle'] = False
          segment.info['patient'] = patient
          segment.info['amplicon'] = amplicon
          segment.info['interval'] = s
          segment.info['segments'] = Set()
          p_interval.append(segment)
        p_interval.sort()
        for cycle in data[patient][1][amplicon]['cycle_map'].keys():
          is_cyclic = data[patient][1][amplicon]['cycle_map'][cycle]['cycle'][0][0] != '0'
          for seg in data[patient][1][amplicon]['cycle_map'][cycle]['cycle']:
            if seg[0:-1] != '0':
              #Find out which intervals it intersects, hopefully only one
              hits = p_interval.intersection([data[patient][1][amplicon]['segment_map'][seg[0:-1]]])
              if len(hits) > 1:
                print "Fucked"
                break    
              elif len(hits) == 0:
                print "Double fucked"
                break
              hits[0][0].info['copy']+=(data[patient][1][amplicon]['cycle_map'][cycle]['copy_count'] *(hits[0][1].end-hits[0][1].start))/ (hits[0][0].end-hits[0][0].start)
              hits[0][0].info['cycle'] = hits[0][0].info['cycle'] or is_cyclic
        intervals.extend(p_interval)
    intervals.sort()          
    for i in intervals:
      hits = [h[0] for h in ensembl.intersection([i])]
      for patient in diseases[disease]:
        #Expression of all transcripts that intersect
        fpkm = sum([data[patient][-1][g.info['gene_id']] for g in hits if g.info['gene_id'] in data[patient][-1]])
        copy_count = [(h[0].end-h[0].start,2**float(h[0].info['seg_mean'])) for h in data[patient][0].intersection([hg19.interval(i.chrom.replace('chr',''), i.start, i.end)])]
        if len(copy_count) > 0:
          copy_count = sum([v[0]*v[1] for v in copy_count])/sum([v[0] for v in copy_count])
        else:
          copy_count = 2
        #Copy count, first check if any intersect
        matches = [h[0] for h in intervals.intersection([i]) if h[0] != i and h[0].info['patient'] == patient]
        #if len(matches) != 0:
        #  copy_count = sum([(v.end-v.start)*v.info['copy'] for v in matches])/sum([(v.end-v.start) for v in matches])
        output.write('%s,%s,%s,%s:%d-%d,%s,%s,%s,%f,%f\n'% (disease,patient,i.info['amplicon'],
                      i.chrom,i.start,i.end,
                     patient==i.info['patient'], patient==i.info['patient'] or len(matches) > 0, i.info['cycle'] if i.info['patient']==patient 
                    else len([m for m in matches if m.info['cycle'] == True]) > 0, fpkm, copy_count))
  output.close()

        
def load_all_amplicons():
  data = dict()
  ks = data_map.keys()
  #ks = [k for k in data_map.keys() if data_map[k][key_map['WGS__cases__project__project_id']] == disease]
  done = 0
  for k in ks:
    res=data_map[k]
    tumor_type = res[key_map['WGS__sampleType']]
    disease_type = res[key_map['WGS__cases__project__project_id']]
    out_dir = '%s/analyses/tumor/%s/min_cnv1/%s'  % (PANCANCER_DIR, disease_type, res[key_map['outdir_basename']])          
    if tumor_type != 'TP':
      continue
    cnv_file = '%s/%s.cnv.bed' % (out_dir,res[key_map['WGS__patient_id']])
    prefix = "%s-MINCNV1" % res[key_map['Tag']]
    summary_file = "%s/output/%s_summary.txt" % (out_dir, prefix)
    fpkm_file = "%s/fpkm_uq.csv" % (out_dir)    
    if res[key_map['aa_finished']] == 'failed':
      #print "Failed %s" % summary_file
      continue
    if not os.path.exists(summary_file) or not os.path.exists(cnv_file) or not os.path.exists(fpkm_file):
      #print "Missing %s" % summary_file
      continue      
    done+=1    
    fpkm = load_fpkm(fpkm_file)          
    cnv = load_bed_file(cnv_file, header=True, map={'chr':1,'start':2,'end':3,'num_probes':4,'seg_mean':5})
    all_map = load_amplicon_file(res)
    targets = load_bed_file("%s/output/target.bed" % out_dir)    
    #Build cycle map list
    hits = [h[1] for h in targets.intersection(cnv)]
    intervals = hg19.interval_list()
    segments = hg19.interval_list()
    cycle_segs = hg19.interval_list()
    cycle_map = {}
    if all_map is not None:
      for amplicon_id in all_map.keys():
        if amplicon_id == 'amplicons':
          continue
        for i in all_map[amplicon_id]['interval_map'].values():
          i.info['amplicon_id'] = amplicon_id
          intervals.append(i)
        for i in all_map[amplicon_id]['segment_map'].values():
          i.info['amplicon_id'] = amplicon_id
          segments.append(i)          
        for (cid, cycles) in all_map[amplicon_id]['cycle_map'].items():
#           if cycles['copy_count'] < 1 or (cycles['cycle'][0][0] == '0' and len(cycles['cycle']) == 3):
#             continue
          ints = [all_map[amplicon_id]['segment_map'][c[0:-1]] for c in cycles['cycle'] if c[0] != '0']
#           if sum([t.end-t.start for t in ints]) < 1000:
#             continue
          foo = [i.info.setdefault('ids',[]).append((amplicon_id,cid)) for i in ints]
          foo = cycle_map.setdefault(amplicon_id,{}).setdefault(cid, hg19.interval_list(ints))
          cycle_map[amplicon_id][cid].sort()
          cycle_segs.extend(ints)
          cycle_segs.sort()
    else:
      continue
    graph = {'intervals':intervals,'segments':segments,'cycles':cycle_segs,'cycle_map':cycle_map}
    data[k] = (cnv, all_map, targets,graph,fpkm)
  #Now for a given amplicon, check to see if they have a path through the breakpoint graph that 
  #satisfies the requirements (i.e., consists with more than an amplified interval or contains a cyclic path
  total_params = 17
  #Now go print out statistics  
  output = open('%s/analyses/aa_overview.csv' % PANCANCER_DIR,'w')  
  output_map = {0:"Sample",1:"Disease",2:"Amplicon",3:"total_amplicon_size",4:"oncogenes", 
                5:"total_amplified_amplicon_size", 6:"average_amplified_copy_count", 
                7:"total_coverage_shifts_with_breakpoint", 8:"total_breakpoint_edges", 
                9:"total_amplicon_chromosome", 10:"total_sequence_edges", 11:"longest_bps_cycle", 
                12:"longest_segments_cycle", 13:"cyclic_cycles",
                14:"total_coverage_shifts", 15:"total_intervals", 16:"amplicon_url"}
  output.write(",".join([output_map[x] for x in xrange(0,total_params)]))
  output.write('\n')
  for k in data.keys():
    res=data_map[k]    
    tumor_type = res[key_map['WGS__sampleType']]
    disease_type = res[key_map['WGS__cases__project__project_id']]
    out_dir = '%s/analyses/tumor/%s/min_cnv1/%s'  % (PANCANCER_DIR, disease_type, res[key_map['outdir_basename']])
    (cnv, all_map, targets,graph) = data[k]
    amplicons = [a for a in data[k][1].keys() if a != 'amplicons']
    for amplicon_id in amplicons:
      oncogenes = len([o for o in data[k][1]['amplicons'][amplicon_id]['OncogenesAmplified'].strip().split(',') if o !=''])
      total_amplicon_size = data[k][1]['amplicons'][amplicon_id]['TotalIntervalSize']
      total_amplified_amplicon_size = data[k][1]['amplicons'][amplicon_id]['AmplifiedIntervalSize']
      average_amplified_copy_count = data[k][1]['amplicons'][amplicon_id]['AverageAmplifiedCopyCount']
      total_coverage_shifts_with_breakpoint = data[k][1]['amplicons'][amplicon_id]['#CoverageShiftsWithBreakpointEdges']
      total_coverage_shifts = data[k][1]['amplicons'][amplicon_id]['#CoverageShifts']
      total_breakpoint_edges = data[k][1]['amplicons'][amplicon_id]['#BreakpointEdges']
      total_amplicon_chromosome = data[k][1]['amplicons'][amplicon_id]['#Chromosomes']
      total_sequence_edges = data[k][1]['amplicons'][amplicon_id]['#SeqenceEdges']      
      total_intervals = data[k][1]['amplicons'][amplicon_id]['#Intervals']     
      #Now cycle specific parameters
      longest_bps_cycle = 'NA'
      longest_segments_cycle = 'NA'
      cyclic_cycles = 'NA'
      prefix = "%s-MINCNV1" % res[key_map['Tag']]
      amplicon_image = "%s/output/%s_amplicon%d.png" % (out_dir,prefix,amplicon_id)
      amplicon_file = "%s/output/%s_amplicon%d_cycles.txt" % (out_dir,prefix,amplicon_id)
      if amplicon_id in graph['cycle_map']:
        cycles = [(cycle_id,sum([c.end-c.start for c in graph['cycle_map'][amplicon_id][cycle_id]]), len([c.end-c.start for c in graph['cycle_map'][amplicon_id][cycle_id]])) for cycle_id in graph['cycle_map'][amplicon_id]]            
        cycles = sorted(cycles, key=lambda x: (x[1],x[2]), reverse=True)
        longest_bps_cycle = cycles[0][1]
        cycles = sorted(cycles, key=lambda x: (x[2],x[1]), reverse=True)
        longest_segments_cycle = cycles[0][2]
        cyclic_cycles = len([cycle_id for cycle_id in graph['cycle_map'][amplicon_id] if data[k][1][amplicon_id]['cycle_map'][cycle_id]['cycle'][0][0] != '0'])
      output_values = {0:k,1:disease_type,2:str(amplicon_id),3:total_amplicon_size,4:str(oncogenes), 
                5:total_amplified_amplicon_size, 6:average_amplified_copy_count, 
                7:total_coverage_shifts_with_breakpoint, 8:total_breakpoint_edges, 
                9:total_amplicon_chromosome, 10:total_sequence_edges, 11:str(longest_bps_cycle), 
                12:str(longest_segments_cycle), 13:str(cyclic_cycles), 14:total_coverage_shifts,
                15:total_intervals, 
                16:"http://genomequery.ucsd.edu:8000/?start=start&cycle=%s&image=%s&prefix=%s" % (amplicon_file, amplicon_image, k)}              
      output.write(",".join([output_values[x] for x in xrange(0,total_params)]))
      output.write('\n')      
  output.close()        
    
    
    
    
  
      
  #Now for a given an amplicon, start filtering 
  Counter([len(data[k][3]['intervals']) for k in data.keys()])
  Counter([len(data[k][1][amp_id]['interval_map'].keys()) for k in data.keys() for amp_id in data[k][1].keys() if amp_id != 'amplicons'])
  len([k for k in data.keys() if len(data[k][3]['cycle_map'].keys()) > 0])
  
  #Filter out cycles that are less than 1
  for k in data.keys():
    (cnv, all_map, targets, graph) = data[k]    
       
  for k in data.keys():
    (cnv, all_map, targets, graph) = data[k]    
      foo = [intervals.append(i) for amp in all_map.keys() if amp != 'amplicons' for i in all_map[amp]['interval_map'].values()]
      foo = [segments.append(i) for amp in all_map.keys() if amp != 'amplicons' for i in all_map[amp]['segment_map'].values()]
      intervals.sort()
      segments.sort()

def determine_output_of_aa():
  data = dict()
  ks = data_map.keys()
  ks = [k for k in data_map.keys() if data_map[k][key_map['WGS__cases__project__project_id']] == disease]
  done = 0
  output = open('%s/analyses/aa_overview.csv' % PANCANCER_DIR,'w')  
  output.write("Sample,Disease,Interval,Size,Genes,Amplification,NumberOfPaths"+
               ",SizeOfLargestAmplicon,IntervalsInLargestAmplicon,SequenceEdgesInLargestAmplicon\n")
  for k in ks:
    res=data_map[k]
    tumor_type = res[key_map['WGS__sampleType']]
    disease_type = res[key_map['WGS__cases__project__project_id']]
    out_dir = '%s/analyses/tumor/%s/min_cnv1/%s'  % (PANCANCER_DIR, disease_type, res[key_map['outdir_basename']])          
    if tumor_type != 'TP':
      continue
    cnv_file = '%s/%s.cnv.bed' % (out_dir,res[key_map['WGS__patient_id']])
    prefix = "%s-MINCNV1" % res[key_map['Tag']]
    summary_file = "%s/output/%s_summary.txt" % (out_dir, prefix)
    if res[key_map['aa_finished']] == 'failed':
      #print "Failed %s" % summary_file
      continue
    if not os.path.exists(summary_file) or not os.path.exists(cnv_file):
      #print "Missing %s" % summary_file
      continue      
    done+=1    
    cnv = load_bed_file(cnv_file, header=True, map={'chr':1,'start':2,'end':3,'num_probes':4,'seg_mean':5})
    all_map = load_amplicon_file(res)
    targets = load_bed_file("%s/output/target.bed" % out_dir)    
    #Build cycle map list
    cycle_segs = hg19.interval_list()
    for amplicon_id in all_map.keys():
      if amplicon_id == 'amplicons':
        continue
      for (cid, cycles) in all_map[amplicon_id]['cycle_map'].items():
        ints = [all_map[amplicon_id]['segment_map'][c[0:-1]] for c in cycles['cycle'] if c[0] != '0']
        ints.sort()
        foo = [i.info.setdefault('ids',[]).append((amplicon_id,cid)) for i in ints]
        cycle_segs.extend(ints)
    cycle_segs.sort()
    hits = [h[1] for h in targets.intersection(cnv)]
    data[k] = (cnv, all_map, targets)
    intervals = hg19.interval_list()
    if all_map is not None:
      foo = [intervals.append(i) for amp in all_map.keys() if amp != 'amplicons' for i in all_map[amp]['interval_map'].values()]
    intervals.sort()
    #Now, for each hit, check to see if it was included in a reconstruction
    for h in hits:
      h_copy = copy.copy(h)
      h_copy.chrom = 'chr%s' % h.chrom
      if all_map is None or len(cycle_segs.intersection([h_copy]))==0:
        output.write("%s,%s,%s,%s,%s,%f,%s,%s,%s,%s\n" % (res[0],disease_type,"%s:%d-%d" % (h.chrom,h.start,h.end), h.end-h.start, genes, 2**float(h.info['seg_mean']), 'NA', 'NA', 'NA', 'NA'))
      else:      
        matches = []  
        foo = [matches.extend(m[0].info['ids']) for m in cycle_segs.intersection([h_copy])]
        matches = sorted(matches, key=lambda x: int(all_map['amplicons'][x[0]]['AmplifiedIntervalSize']), reverse=True)
        genes = len(ensembl.intersection([h_copy]))
      output.write("%s,%s,%s,%s,%s,%f,%s,%s,%s,%s\n" % (res[0],disease_type,"%s:%d-%d" % (h.chrom,h.start,h.end), h.end-h.start, genes, 2**float(h.info['seg_mean']), len(matches), all_map['amplicons'][matches[0][0]]['AmplifiedIntervalSize'], all_map['amplicons'][matches[0][0]]['#Intervals'], all_map['amplicons'][matches[0][0]]['#SeqenceEdges']))
  output.close()
        
    
    #For each 
    
def get_random_segment(chr=None, length=1, replicates=1):
  intervals = hg19.interval_list([])
  if chr is None:
    chrLenTotal = sum(hg19.chrLen.values())
    while replicates > 0:
      replicates = replicates - 1      
      (chr,pos) = hg19.chrPos(random.randint(1,chrLenTotal))
      intervals.append(hg19.interval(chr,max(1,pos-int(length/2)),pos+int(length/2+0.5)))
  return intervals    

def load_fish_file(file='%s/analyses/FISHfile.txt' % (PANCANCER_DIR)):
  data = {}
  current_sample = None
  for line in open(file, 'r'):
    res = line.strip().split('\t')
    if res[1] == '':
      res[1] = current_sample
    else:
      current_sample = res[1]
    data.setdefault(current_sample,{})[res[0]] = res
  return data

def graph_statistics():
  data = load_fish_file()
  read_counts = dict()
  ks = keeper.keys()
  #ks = [k for k in complete.keys() if complete[k][key_map['WGS__cases__project__project_id']] == disease]
  for k in ks:
    res=keeper[k]
    tumor_type = res[key_map['WGS__sampleType']]
    disease_type = res[key_map['WGS__cases__project__project_id']]
    out_dir = '%s/analyses/tumor/%s/min_cnv1/%s/'  % (PANCANCER_DIR, disease_type, res[key_map['outdir_basename']])          
    if tumor_type != 'TP':
      continue
    cnv_bed_file =  '%s/%s.hg19.bed' % (out_dir,res[key_map['WGS__patient_id']])
    if not os.path.exists(cnv_bed_file):
      continue
    cnv = load_cnv_file(cnv_bed_file)
    all_map = load_amplicon_file(res)
    read_counts[k] = (cnv,all_map)
  #Compute FPKM for any segment of a cycle that is cyclic
  output = open('%s/analyses/graph.csv' % (PANCANCER_DIR),'w')
  output.write('Sample,Disease,Amplicon_id,total_copy_count,total_segment_count,total_length'+
        ',total_chrom,total_paths,total_circular_cycles,total_breakpoints,'+
        'longest_path,size_of_longest_path,longest_heaviest_path,size_of_heaviest_path,oncogenes,type\n')  
  for k in complete.keys():
    res = complete[k]
    if k not in read_counts:
      continue
    #Compute for each amplicon
    for (amplicon_id, all_map) in read_counts[k][1].items():
      if amplicon_id == 'amplicons':
        continue
      cycles = [(c,cycle) for (c,cycle) in all_map['cycle_map'].items()]      
      if len(cycles) == 0:
        continue      
      cycles = sorted(cycles, key=lambda x: x[1]['copy_count'], reverse=True)      
      total_cycle_count = sum([c['copy_count'] for c in all_map['cycle_map'].values()])      
      total_segment_count = len(all_map['segment_map'].keys())
      total_length = sum([seg.end-seg.start for seg in all_map['segment_map'].values()])
      total_chrom = len(Set([seg.chrom for seg in all_map['segment_map'].values()]))
      total_cycles = len(all_map['cycle_map'].values())
      total_circular_cycles =  len([v for v in all_map['cycle_map'].values() if v['cycle'][0][0] != '0'])
      total_breakpoints = sum([len([x for x in v['cycle'] if x[0] != '0'])-1 for v in all_map['cycle_map'].values()])
      longest_dominate_path = len([c for c in cycles[0][1]['cycle'] if c[0] != '0'])
      size_of_dominate_path = sum([all_map['segment_map'][c[0:-1]].end-all_map['segment_map'][c[0:-1]].start for c in cycles[0][1]['cycle'] if c[0] != '0'])      
      cycles = sorted(cycles, key=lambda x: (len([c for c in x[1]['cycle'] if c[0] != '0']), sum([all_map['segment_map'][c[0:-1]].end-all_map['segment_map'][c[0:-1]].start for c in x[1]['cycle'] if c[0] != '0'])), reverse=True)      
      longest_path = len([c for c in cycles[0][1]['cycle'] if c[0] != '0'])
      size_of_longest_path = sum([all_map['segment_map'][c[0:-1]].end-all_map['segment_map'][c[0:-1]].start for c in cycles[0][1]['cycle'] if c[0] != '0'])            
      oncogenes = len([o for o in read_counts[k][-1]['amplicons'][amplicon_id]['OncogenesAmplified'].split(',') if o != ''])
      output.write('%s,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s\n' % (res[0], res[5], amplicon_id, total_cycle_count,total_segment_count,total_length,total_chrom,total_cycles,total_circular_cycles,total_breakpoints,longest_path,size_of_longest_path,longest_dominate_path,size_of_dominate_path,oncogenes,'TCGA'))
  #Now load up Viraj's results
  turner = load_sample_list()
  turner_dir = '/pedigree/projects/extrachromosome/data/turner2017/reconstruction/run14'
  for (key, value) in turner.items():
    all_map = load_amplicon_file(res=None, summary_file='%s/%s_summary.txt' % (turner_dir,key), prefix=key, out_dir = turner_dir)
    value.append(all_map)
  
  for k in turner.keys():
    #Compute for each amplicon
    for (amplicon_id, all_map) in turner[k][-1].items():
      if amplicon_id == 'amplicons':
        continue    
      cycles = [(c,cycle) for (c,cycle) in all_map['cycle_map'].items()]      
      if len(cycles) == 0:
        continue      
      cycles = sorted(cycles, key=lambda x: x[1]['copy_count'], reverse=True)      
      total_cycle_count = sum([c['copy_count'] for c in all_map['cycle_map'].values()])      
      total_segment_count = len(all_map['segment_map'].keys())
      total_length = sum([seg.end-seg.start for seg in all_map['segment_map'].values()])
      total_chrom = len(Set([seg.chrom for seg in all_map['segment_map'].values()]))
      total_cycles = len(all_map['cycle_map'].values())
      total_circular_cycles =  len([v for v in all_map['cycle_map'].values() if v['cycle'][0][0] != '0'])
      total_breakpoints = sum([len([x for x in v['cycle'] if x[0] != '0'])-1 for v in all_map['cycle_map'].values()])
      longest_dominate_path = len([c for c in cycles[0][1]['cycle'] if c[0] != '0'])
      size_of_dominate_path = sum([all_map['segment_map'][c[0:-1]].end-all_map['segment_map'][c[0:-1]].start for c in cycles[0][1]['cycle'] if c[0] != '0'])      
      cycles = sorted(cycles, key=lambda x: (len([c for c in x[1]['cycle'] if c[0] != '0']), sum([all_map['segment_map'][c[0:-1]].end-all_map['segment_map'][c[0:-1]].start for c in x[1]['cycle'] if c[0] != '0'])), reverse=True)      
      longest_path = len([c for c in cycles[0][1]['cycle'] if c[0] != '0'])
      size_of_longest_path = sum([all_map['segment_map'][c[0:-1]].end-all_map['segment_map'][c[0:-1]].start for c in cycles[0][1]['cycle'] if c[0] != '0'])            
      oncogenes = len([o for o in turner[k][-1]['amplicons'][amplicon_id]['OncogenesAmplified'].split(',') if o != ''])
      type = 'Unknown' if turner[k][1] not in data else data[turner[k][1]]
      if type != 'Unknown':
        type = [v[6].split()[0].split(':')[-1] for v in data[turner[k][1]].values() if len(v) > 6]
        if len(type) == 0:
          type = 'Unknown'
#         elif len(type) > 1:
#           print 'Fuck'
#           print type
        else:
          type = type[0]
      output.write('%s,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s\n' % (k, turner[k][3].strip(), amplicon_id, total_cycle_count,total_segment_count,total_length,total_chrom,total_cycles,total_circular_cycles,total_breakpoints,longest_path,size_of_longest_path,longest_dominate_path,size_of_dominate_path,oncogenes,type))
  output.close()
  

def load_sample_list(sample_list = '/pedigree/projects/extrachromosome/data/turner2017/sample_list.csv'):
  input = open(sample_list, 'r')
  sample_map = {}
  for line in input:
    res = line.strip().split(',')
    sample_map[res[0]] = res
  return sample_map

def simulate_sample():
  read_counts = dict()
  ks = complete.keys()
  #ks = [k for k in complete.keys() if complete[k][key_map['WGS__cases__project__project_id']] == disease]
  for k in ks:
    res=complete[k]
    tumor_type = res[key_map['WGS__sampleType']]
    disease_type = res[key_map['WGS__cases__project__project_id']]
    out_dir = '%s/analyses/tumor/%s/min_cnv1/%s/'  % (PANCANCER_DIR, disease_type, res[key_map['outdir_basename']])          
    if tumor_type != 'TP':
      continue
    flagstat_file = "%s/flagstat.csv" % (out_dir)
    read_file = "%s/read_counts.csv" % (out_dir)
    fpkm_file = "%s/fpkm_uq.csv" % (out_dir)
    if not os.path.exists(flagstat_file) or not os.path.exists(read_file) or not os.path.exists(fpkm_file):
      continue
    flagstat = read_samtools_flagstat(flagstat_file)
    reads = read_readcounts(read_file)
    read_counts[k] = (reads, flagstat)            
    cnv_bed_file =  '%s/%s.hg19.bed' % (out_dir,res[key_map['WGS__patient_id']])
    cnv = load_cnv_file(cnv_bed_file)
    fpkm = load_fpkm(fpkm_file)
    all_map = load_amplicon_file(res)
    read_counts[k] = (reads, flagstat, cnv, fpkm, all_map)
  #Compute FPKM for any segment of a cycle that is cyclic
  output = open('%s/analyses/graph.csv' % (PANCANCER_DIR),'w')
  output.write('Sample,Disease,Amplicon_id,total_cycle_count,total_segment_count,total_length'+
        ',total_chrom,total_cycles,total_circular_cycles,total_breakpoints,total_intra_edges,'+
        'longest_path,size_of_longest_path,longest_dominate_path,size_of_dominate_path\n')  
  for k in complete.keys():
    res = complete[k]
    if k not in read_counts:
      continue
    #Compute for each amplicon
    for (amplicon_id, all_map) in read_counts[k][4].items():
      cycles = [(c,cycle) for (c,cycle) in all_map['cycle_map'].items()]      
      if len(cycles) == 0:
        continue      
      cycles = sorted(cycles, key=lambda x: x[1]['copy_count'], reverse=True)      
      total_cycle_count = sum([c['copy_count'] for c in all_map['cycle_map'].values()])      
      total_segment_count = len(all_map['segment_map'].keys())
      total_length = sum([seg.end-seg.start for seg in all_map['segment_map'].values()])
      total_chrom = len(Set([seg.chrom for seg in all_map['segment_map'].values()]))
      total_cycles = len(all_map['cycle_map'].values())
      total_circular_cycles =  len([v for v in all_map['cycle_map'].values() if v['cycle'][0][0] != '0'])
      total_breakpoints = sum([len([x for x in v['cycle'] if x[0] != '0'])-1 for v in all_map['cycle_map'].values()])
      total_intra_edges = 0      
      for (c,cycle) in all_map['cycle_map'].items():
        current_chrom = None
        for seg in cycle['cycle']:
          if current_chrom == None:
            current_chrom = segment.chrom
          else:
            if current_chrom != segment.chrom:
              total_intra_edges+=1
            current_chrom = segment.chrom
      longest_dominate_path = len([c for c in cycles[0][1]['cycle'] if c[0] != '0'])
      size_of_dominate_path = sum([all_map['segment_map'][c[0:-1]].end-all_map['segment_map'][c[0:-1]].start for c in cycles[0][1]['cycle'] if c[0] != '0'])      
      cycles = sorted(cycles, key=lambda x: (len([c for c in x[1]['cycle'] if c[0] != '0']), sum([all_map['segment_map'][c[0:-1]].end-all_map['segment_map'][c[0:-1]].start for c in x[1]['cycle'] if c[0] != '0'])), reverse=True)      
      longest_path = len([c for c in cycles[0][1]['cycle'] if c[0] != '0'])
      size_of_longest_path = sum([all_map['segment_map'][c[0:-1]].end-all_map['segment_map'][c[0:-1]].start for c in cycles[0][1]['cycle'] if c[0] != '0'])            
      output.write('%s,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n' % (res[0], res[5], amplicon_id, total_cycle_count,total_segment_count,total_length,total_chrom,total_cycles,total_circular_cycles,total_breakpoints,total_intra_edges,longest_path,size_of_longest_path,longest_dominate_path,size_of_dominate_path))  
#       number_of_segments = len(all_map['segment_map'].keys())
#       number_of_intra_edges = 0
#       longest_path = 0
#       cycles = [(c,cycle) for (c,cycle) in all_map['cycle_map'].items()]            
#       cycles = sorted(cycles, key=lambda x: x[1]['copy_count'], reverse=True)      
#       total_cycle = sum([c['copy_count'] for c in all_map['cycle_map'].values()])      
#       cycle_explain_data = 0
#       curr_data = 0      
#       for c in cycles:
#         cycle_explain_data+=1
#         curr_data+=c[1]['copy_count']
#         if curr_data/total_cycle >= .80:
#           break          
#       longest_path_of_dominate_cycle = len([c for c in cycles[0][1]['cycle'] if c[0] != '0'])
#       size_of_dominate_cycle = sum([all_map['segment_map'][c[0:-1]].end - all_map['segment_map'][c[0:-1]].start for c in cycles[0][1]['cycle'] if c[0] != '0'])
#       longest_cycle = 0
#       largest_cycle = 0      
#       for (c,cycle) in all_map['cycle_map'].items():
#         is_cycle = True
#         length_of_cycle = 0
#         size_of_cycle = 0
#         current_chrom = None
#         for seg in cycle['cycle']:
#           if seg[0] == '0':
#             is_cycle = False
#           else:
#             length_of_cycle+=1
#             segment = all_map['segment_map'][seg[0:-1]]
#             size_of_cycle+=segment.end-segment.start
#             if current_chrom == None:
#               current_chrom = segment.chrom
#             else:
#               if current_chrom != segment.chrom:
#                 number_of_intra_edges+=1
#               current_chrom = segment.chrom
#         if longest_cycle < length_of_cycle:
#           longest_cycle = length_of_cycle                             
#         if largest_cycle < size_of_cycle:
#           largest_cycle = size_of_cycle
  #Compute FPKM for any segment of a cycle that is cyclic
  output = open('%s/analyses/simulation.csv' % (PANCANCER_DIR),'w')
  #output.write('Sample,Disease,Gene_name ,Gene,Structure,FPKM,CNV,Oncogene\n')  
  counter = -1
  for disease in diseases:
    counter+=1
    #if counter % 6 != index:
    #  continue
    cycle_maps = {}
    segment_maps = {}
    interval_maps = {}
    samples = [v[key_map['Tag']] for v in complete.values() if v[key_map['WGS__cases__project__project_id']] == disease and v[key_map['Tag']] in read_counts]    
    for reference_sample in samples:
      for (reference_amp_id,all_map) in read_counts[reference_sample][4].items():
        cycles = [(cid,cmap) for (cid,cmap) in all_map['cycle_map'].items() if cmap['cycle'][0][0] != '0' and cmap['cycle'][-1][0] != '0']
        interval_maps.setdefault(reference_sample, hg19.interval_list()).extend(all_map['interval_map'].values())
        interval_maps[reference_sample].sort()
        if len(cycles) == 0:
          continue
        for c in cycles:
          cycle_interval_list = hg19.interval_list()        
          for segment in c[1]['cycle']:
            if segment[:-1] == '0':
              continue
            c_seg = all_map['segment_map'][segment[:-1]]
            c_seg.info['cycle_id'] = c[0]
            c_seg.info['cycle_map'] = c[1]            
            cycle_interval_list.append(c_seg)
          cycle_interval_list.sort()
#         if sum([(c.end-c.start) for c in cycle_interval_list]) > 50000:
#           continue
          cycle_maps.setdefault(reference_sample,{}).setdefault(reference_amp_id,{}).setdefault(
                                c[0], hg19.interval_list()).extend(cycle_interval_list)
          segment_maps.setdefault(reference_sample,hg19.interval_list()).extend(cycle_interval_list)
          segment_maps[reference_sample].sort()
    #Find all intersecting genes
    all_intervals = hg19.interval_list([i for v in interval_maps.values() for i in v])
    all_intervals.sort()
    rep_array = [get_random_segment(length=(i.end-i.start), replicates=100) for i in all_intervals]
    for i in xrange(0,len(rep_array[0])):
      array = hg19.interval_list([rep_array[j][i] for j in xrange(0,len(rep_array))])
      for a in array:
        a.chrom='chr%s' % a.chrom      
      array.sort()
      oncogene_count = 0
      gene_count = len([h for h in ensembl.intersection(array)])
      oncogene_count = len([h[0] for h in hg19.oncogene_list.intersection(array)])
      output.write('%s,%d,%d,%d\n' % (disease,i,gene_count,oncogene_count))        
  output.close()  

def load_fusion_data(file="%s/analyses/pancanfus.txt" % (PANCANCER_DIR)):
  input = open(file, 'r')
  header = input.next().strip().split('\t')
  fusion_map = {}  
  for line in input:
    res = line.strip().split('\t')
    disease = "TCGA-%s" % res[0]
    sample = "TCGA-%s-%s" % (res[1].split('.')[0], res[1].split('.')[1])    
    if res[1].split('.')[2] != '01A':
      continue
    intervalA = hg19.interval("%s" % res[4],res[17],res[18]) 
    intervalB = hg19.interval("%s" % res[5],res[19],res[20])     
    fusion_map.setdefault(sample,hg19.interval_list()).extend([intervalA,intervalB])
  for value in fusion_map.values():
    value.sort()
  return fusion_map
  
def check_fusion():
  amplicon_list = {}
  output = open('%s/analyses/fusion.csv' % PANCANCER_DIR,'w')
  output.write('Sample,Disease,Fusion\n')
  for disease in [d for d in os.listdir('%s/analyses/tumor' % PANCANCER_DIR) if d.find('TCGA') != -1]:
    amplicon_list[disease] = collect_amplicons(complete,disease=disease)
    sample_map = {}
    foo = [sample_map.setdefault(a.info['meta'][key_map['WGS__patient_id']], hg19.interval_list()).append(a) for a in amplicon_list[disease]]
    for value in sample_map.values():
      value.sort()
    for key in sample_map.keys():
      hits = sample_map[key].intersection(fusion_map[key] if key in fusion_map else [])
      output.write("%s,%s,%s\n" % (key, disease, len(hits) > 0))
  output.close()
  
def download_fpkm_values(complete, key_map):
  for (k,res) in complete.items():
    tumor_type = res[key_map['WGS__sampleType']]
    disease_type = res[key_map['WGS__cases__project__project_id']]
    out_dir = 'analyses/tumor/%s/min_cnv1/%s'  % (disease_type, res[key_map['outdir_basename']])          
    if tumor_type != 'TP':
      continue
    fpkm_file = res[key_map['rnaseq_fpkm_uq']]
    if os.path.exists("%s/%s/fpkm_uq.csv" % (PANCANCER_DIR, out_dir)):
      continue
    path = tempfile.mkdtemp(dir=TEMP_DIR)
    os.system("$WORK/bin/gdc-client download --http-chunk-size 20971520 -n 1 -d %s %s" % (path, fpkm_file))
    file = [f for f in os.listdir('%s/%s' % (path,fpkm_file)) if f.find('FPKM-UQ') != -1][0]
    os.chdir("%s/%s" % (path,fpkm_file))
    os.system('gzip -d -c %s > fpkm_uq.csv' % file)
    os.system('mv fpkm_uq.csv %s/%s/' % (PANCANCER_DIR, out_dir))
    os.chdir(PANCANCER_DIR)
    os.system('rm -rf %s' % path)

def load_meta_file():
  data_map = dict()
  for disease in [d for d in os.listdir('%s/analyses/tumor' % PANCANCER_DIR) if d.find('TCGA') != -1]:
  #for disease in ['TCGA-GBM']:
  #for disease in ['TCGA-BRCA']:
    file="%s/analyses/tumor/%s/min_cnv1/file_map.txt" % (PANCANCER_DIR, disease)
    input = open(file, 'r')
    header=input.next().split('\t')
    idx = -1
    key_map = dict([(header[i].replace('"',""),i) for i in xrange(0,len(header))])
    i+=1
    key_map['rnaseq_gsc'] = i
    key_map['sample'] = i+1
    key_map['meta'] = i+2  
    key_map['rnaseq_fpkm_uq'] = i+3    
    for i in input:
      res = [r.replace('"','') for r in i.strip().split('\t')]
      data_map.setdefault(res[0],res)
    #Find if the file has an RNA-seq file
  return (data_map, key_map)

def build_amplicon_list(complete):
  amplicon_list = {}
  for disease in [d for d in os.listdir('%s/analyses/tumor' % PANCANCER_DIR) if d.find('TCGA') != -1]:
    amplicon_list[disease] = collect_amplicons(complete,disease=disease)
    run_read_count(complete,disease=disease)


def run_read_count(complete, disease='TCGA-BRCA'):
  output = open('%s/scripts/counter.%s.sh' % (PANCANCER_DIR,disease), 'w')
  keys = complete.keys()
  random.shuffle(keys)
  for k in keys:
    res=complete[k]
    tumor_type = res[key_map['WGS__sampleType']]
    disease_type = res[key_map['WGS__cases__project__project_id']]
    out_dir = 'analyses/tumor/%s/min_cnv1/%s'  % (disease_type, res[key_map['outdir_basename']])          
    if tumor_type != 'TP' or disease_type != disease:
      continue
    rna_bam_file = res[key_map['rnaseq_gsc']]
    try:
      res = subprocess.check_output(['gsutil','ls', '-lh', 'gs://aa-data-repo/%s/read_counts.csv' % (out_dir)],stderr=subprocess.STDOUT)      
      if len(res.strip()) != 0 and float(res.strip().split(' ')[0]) > 0:
        continue
    except subprocess.CalledProcessError as error:
      print "Try to analyze %s" % out_dir
    res = subprocess.check_output(['gsutil','ls', '-lh', rna_bam_file],stderr=subprocess.STDOUT).strip()
    mem = int(float(res.split(' ')[0])+5)
    output.write("dsub --image us.gcr.io/aa-test-175718/read_count --preemptible --project aa-test-175718 --zones \"us-west1-*\" --logging gs://aa-data-repo/%s/logging/ \\\n\t--output LOG_FILE=gs://aa-data-repo/%s/log.txt \\\n\t \\\n\t--output-recursive OUTPUT_DIR=gs://aa-data-repo/%s --disk-size %d \\\n\t--input SEGMENT_FILE=gs://aa-data-repo/metadata/amplicons.%s.csv BAM_FILE=%s BAI_FILE=%s.bai \\\n\t --command 'sh /home/read_count.sh'\n" %  (out_dir,out_dir,out_dir, mem, disease, rna_bam_file,rna_bam_file))
#     path = tempfile.mkdtemp(dir=TEMP_DIR) 
#     os.system('gsutil cp %s.bai %s/rna.bam.bai' % (rna_bam_file, path))    
#     os.system('gsutil cp %s %s/rna.bam' % (rna_bam_file, path))
#     os.system("python %s/python/read_count.py -b %s/rna.bam -s %s/metadata/amplicons.%s.csv -f %s/flagstat.csv -o %s/read_counts.csv"  % (PANCANCER_DIR, path, PANCANCER_DIR, disease, path, path))
#     os.system("samtools flagstat %s/rna.bam > %s/flagstat.csv"  % (path, path))    
#     os.system('gsutil cp %s/*.csv gs://aa-data-repo/%s' % (path, out_dir))
#     os.system('rm -rf %s' % path)    
  output.close()  

def get_rnaseq_bucket(data_map, service, key_map):
  keeper = dict()
  for d in data_map.values():
    try:
      res = service.samples().get(sample_barcode=d[key_map['WGS__portion_id']][0:-3]).execute()
      cases = service.cases().get(case_barcode=d[key_map['WGS__portion_id']][0:-7]).execute()    
    except:
      continue
    wgs = None
    rnaseq = None
    if 'data_details' not in res:
      continue
    for r in res['data_details']:   
      if 'experimental_strategy' not in r:
        continue
      if r['experimental_strategy'] == 'WGS' and wgs is None:
        wgs = r
      if r['experimental_strategy'] == 'RNA-Seq' and rnaseq is None:
        rnaseq = r
    if wgs is not None and rnaseq is not None and wgs['file_name_key'] == d[key_map['WGS__file_gcs_url']]:
      tmp = copy.copy(d)
      tmp.append(rnaseq['file_name_key'])
      tmp.append(res)      
      tmp.append(cases)            
      tmp.append(None)
      keeper[d[0]] = tmp
    else:
      tmp = copy.copy(d)
      tmp.append(None)
      tmp.append(res)      
      tmp.append(cases)    
      tmp.append(None)        
      keeper[d[0]] = tmp
  idx = 0
  complete = copy.copy(keeper)
  for k in keeper.keys():
    d = keeper[k]
    if len(complete[k]) < 80:
      del complete[k]
      del keeper[k]
      continue      
    if not os.path.exists('%s/analyses/tumor/%s/min_cnv1/%s/output/finish_flag.txt'  % (PANCANCER_DIR, d[key_map['WGS__cases__project__project_id']], d[key_map['outdir_basename']])) or d[key_map['aa_finished']] != "success":
      del complete[k]
      del keeper[k]
      continue
    if complete[k][key_map['rnaseq_gsc']] is None:
      del complete[k]
    print "Working on %s" % d[key_map['outdir_basename']]
    idx+=1
  print idx  
  return(complete,keeper)  

def fix_cnv_file(cnv_file,output_file):
  os.system("awk '{print \"chr\"$2,$3,$4,$5,$6}' %s | tail -n +2 > %s.tmp" % (cnv_file, cnv_file))
  os.system('$WORK/programs/liftOver/liftOver %s.tmp $WORK/programs/liftOver/hg38ToHg19.over.chain %s %s.unmapped' % (cnv_file, output_file, output_file))
  os.system('rm %s.tmp' % (cnv_file))

def compute_flagstat():
  for k in complete.keys():
    res=complete[k]
    tumor_type = res[key_map['WGS__sampleType']]
    disease_type = res[key_map['WGS__cases__project__project_id']]
    out_dir = '%s/analyses/tumor/%s/min_cnv1/%s/'  % (PANCANCER_DIR, disease_type, res[key_map['outdir_basename']])          
    if tumor_type != 'TP':
      continue
    rna_bam_file = res[key_map['rnaseq_gsc']]      
    for chr in chrs:
      os.system('samtools view -hb %s|samtools flagstat - > %s/%s.flagstat 2>&1' % (rna_bam_file, out_dir, chr))

def read_cycle_file(cycle_file):
  input = open(cycle_file)
  segment_map = {}
  interval_map = {}
  cycle_map = {}
  for line in input:
    res = line.strip().split('\t')
    if res[0] == 'Interval':
      interval_map[res[1]] = hg19.interval(res[2] if res[2].find('chr') != -1 else "chr%s" % res[2],int(res[3]),int(res[4]),info={'line':res})
    if res[0] == 'Segment':
      segment_map[res[1]] = hg19.interval(res[2] if res[2].find('chr') != -1 else "chr%s" % res[2],int(res[3]),int(res[4]),info={'line':res})
    if res[0].find('Cycle') != -1:
      segments = res[0].split(';')
      cycle_name = segments[0].split('=')[-1]
      copy_count = float(segments[1].split('=')[-1])
      cycle = segments[2].split('=')[-1].split(',')
      if cycle[0][0] != '0' and len([c for c in cycle if c[0] == '0']) > 0:
        idx = cycle.index([c for c in cycle if c == '0+'][0])
        cycle = cycle[idx:] + cycle[:idx]
      cycle_map[cycle_name] = {'cycle':cycle,'copy_count':copy_count}
  input.close()
  return (segment_map, interval_map, cycle_map)    

def load_amplicon_file(res, summary_file = None, prefix = "", out_dir = ""):
  if summary_file == None:
    disease_type = res[key_map['WGS__cases__project__project_id']]
    out_dir = '%s/analyses/tumor/%s/min_cnv1/%s/output'  % (PANCANCER_DIR, disease_type, res[key_map['outdir_basename']])              
    prefix = "%s-MINCNV1" % res[key_map['Tag']]
    summary_file = "%s/%s_summary.txt" % (out_dir, prefix)
  all_map = {}
  if os.path.exists(summary_file):    
    amps = read_summary_file(summary_file)
    if amps is None:
      return None
  for (ai, a) in amps.items():
    ints = a['Intervals'].split(',')
    cycle_file = "%s/%s_amplicon%s_cycles.txt" % (out_dir,prefix,ai)
    (segment_map, interval_map, cycle_map) = read_cycle_file(cycle_file)
    all_map.setdefault(ai,{})['segment_map'] = segment_map
    all_map.setdefault(ai,{})['interval_map'] = interval_map
    all_map.setdefault(ai,{})['cycle_map'] = cycle_map
  all_map['amplicons'] = amps
  return all_map

def collect_amplicons(complete,disease='TCGA-GBM'):
  amplicons = hg19.interval_list()
  for k in complete.keys():
    res=complete[k]
    tumor_type = res[key_map['WGS__sampleType']]
    disease_type = res[key_map['WGS__cases__project__project_id']]
    out_dir = '%s/analyses/tumor/%s/min_cnv1/%s/'  % (PANCANCER_DIR, disease_type, res[key_map['outdir_basename']])              
    if tumor_type != 'TP' or disease_type != disease:
      continue
    cnv_bed_file =  '%s/%s.hg19.bed' % (out_dir,res[key_map['WGS__patient_id']])
    cnv = load_cnv_file(cnv_bed_file)      
    summary_file = "%s/output/%s-MINCNV1_summary.txt" % (out_dir,res[key_map['Tag']])    
    if os.path.exists(summary_file):    
      amps = read_summary_file(summary_file)
    #Now for each amplicon add meta data on its segments and cycles
    for (ai, a) in amps.items():
      ints = a['Intervals'].split(',')
      cycle_file = "%s/output/%s-MINCNV1_amplicon%s_cycles.txt" % (out_dir,res[key_map['Tag']],ai)
      (segment_map, interval_map, cycle_map) = read_cycle_file(cycle_file)
      a['segment_map'] = segment_map
      a['interval_map'] = interval_map
      a['cycle_map'] = cycle_map
      for i in ints:     
        it = hg19.interval(i,info={'amplicon':a, 'amplicon_id':ai, 'meta':res})
        amplicons.append(it)
  amplicons.sort()
  amplicons = hg19.interval_list([a for a in amplicons if a.chrom in chrs])
#   output = open('%s/metadata/amplicons.%s.csv' % (PANCANCER_DIR,disease),'w')
#   output.write('Sample,Disease,Amplicon_id,Segment\n')
#   for a in amplicons:
#     id = a.info['meta'][0]
#     disease = a.info['meta'][key_map['WGS__cases__project__project_id']]
#     amp = a.info['amplicon_id']
#     seg = "chr%s:%d-%d" % (a.chrom, a.start, a.end)
#     output.write('%s,%s,%s,%s\n' % (id, disease, amp, seg))
#   output.close()
#   
#   output = open('%s/metadata/complete.%s.csv' % (PANCANCER_DIR,disease),'w')
#   output.write('Sample,Disease,Bucket\n')
#   for k in complete.keys():
#     res=complete[k]
#     id = res[0]
#     tumor_type = res[key_map['WGS__sampleType']]
#     disease_type = res[key_map['WGS__cases__project__project_id']]
#     rna_bam_file = res[key_map['rnaseq_gsc']]
#     if tumor_type != 'TP':
#       continue    
#     output.write('%s,%s,%s\n' % (id, disease, rna_bam_file))
#   output.close()
  return amplicons


def read_readcounts(read_file):
  read_count_intervals = {}
  input = open(read_file, 'r')
  for line in input:
    res = line.strip().split(',')  
    foo = read_count_intervals.setdefault("%s_%s_%s" % (res[0],res[2],res[1] if res[1].find('chr') != -1 else "chr%s" % res[1]),res)
  return read_count_intervals

def read_samtools_flagstat(flagstat_file):
  flagstat = {}
  input = open(flagstat_file,'r')  
  res = input.next().strip().split(' ')
  flagstat['total_reads'] = (int(res[0]),int(res[2]))
  res = input.next().strip().split(' ')
  res = input.next().strip().split(' ')
  flagstat['total_mapped'] = (int(res[0]),int(res[2]))
  res = input.next().strip().split(' ')
  flagstat['total_paired_sequenced'] = (int(res[0]),int(res[2]))
  res = input.next().strip().split(' ')
  flagstat['read1'] = (int(res[0]),int(res[2]))
  res = input.next().strip().split(' ')
  flagstat['read2'] = (int(res[0]),int(res[2]))
  res = input.next().strip().split(' ')
  flagstat['properly_paired'] = (int(res[0]),int(res[2]))
  res = input.next().strip().split(' ')
  flagstat['with_itself_and_mate'] = (int(res[0]),int(res[2]))
  res = input.next().strip().split(' ')
  flagstat['singletons'] = (int(res[0]),int(res[2]))
  res = input.next().strip().split(' ')
  flagstat['mate_to_diff_chr'] = (int(res[0]),int(res[2]))
  res = input.next().strip().split(' ')
  flagstat['mate_to_diff_chr_high_q'] = (int(res[0]),int(res[2]))
  return flagstat    

def parse_segment_file(segment_file):
  input = open(segment_file, 'r')
  amplicons = hg19.interval_list()
  line = input.next()
  for line in input:
    res = line.strip().split(',')
    amplicons.append(hg19.interval(res[3],info={'sample':res[0],'disease':res[1],'id':res[2]}))
  return amplicons

def samtools_flagstat(bam_file, output_file, threads=4):
  os.system('samtools flagstat -@ %d  %s > %s 2>&1'% (threads, bam_file, output_file))
# 
#   data = {}
#   for k in complete.keys():
#     res=complete[k]
#     tumor_type = res[key_map['WGS__sampleType']]
#     disease_type = res[key_map['WGS__cases__project__project_id']]
#     out_dir = '%s/analyses/tumor/%s/min_cnv1/%s/'  % (PANCANCER_DIR, disease_type, res[key_map['outdir_basename']])          
#     if tumor_type != 'TP':
#       continue
#     os.system('')
#     rna_bam_file = res[key_map['rnaseq_gsc']]
#     rna_bam = pysam.Samfile(rna_bam_file, 'rb')     
    

def compute_egf_mortality(keeper, complete, disease='TCGA-BRCA'):
  egfr = hg19.interval_list([g for g in hg19.gene_list if g.info['Name']=='EGFR'])
  output = open('%s/analyses/egf_survival.%s.csv' % (PANCANCER_DIR, disease),'w')
  output.write('Sample,Disease,EGFR,Death,Total_Size,Breakpoint_Edges,Oncogene,Status,Relapse\n')
  for k in keeper.keys():  
    res=keeper[k]
    tumor_type = res[key_map['WGS__sampleType']]
    disease_type = res[key_map['WGS__cases__project__project_id']]
    relapse = "NA" if 'new_tumor_event_after_initial_treatment' not in res[key_map['meta']]['clinical_data'] else res[key_map['meta']]['clinical_data']['new_tumor_event_after_initial_treatment']
    if tumor_type != 'TP':
      continue
    out_dir = '%s/analyses/tumor/%s/min_cnv1/%s/'  % (PANCANCER_DIR, disease_type, res[key_map['outdir_basename']])      
    cnv_file = '%s/%s.cnv' % (out_dir,res[key_map['WGS__patient_id']])    
    cnv_bed_file =  '%s/%s.hg19.bed' % (out_dir,res[key_map['WGS__patient_id']])
    if not os.path.exists(cnv_file):
      os.system('gsutil cp %s %s' % (res[key_map['CNV__file_gcs_url']],cnv_file))
    if not os.path.exists(cnv_bed_file):
      fix_cnv_file(cnv_file, cnv_bed_file)
    cnv = load_cnv_file(cnv_bed_file)
    hits = Set([h[1] for h in egfr.intersection(cnv)])
    if len(hits) == 0:
      amp = 1
    else:
      amp = sum([h.info['Amplification'] for h in hits])/len(hits)
    summary_file = "%s/output/%s-MINCNV1_summary.txt" % (out_dir,res[key_map['Tag']])    
    if os.path.exists(summary_file):    
      amplicons = read_summary_file(summary_file)
      output.write("%s,%s,%f,%s,%d,%d,%s,%s,%s\n" % (res[key_map['WGS__patient_id']],
                                disease_type,amp,  
                                "%d" % res[key_map['meta']]['clinical_data']['days_to_last_known_alive'] if 'days_to_last_known_alive' in res[key_map['meta']]['clinical_data'] else 'NA',
                                sum([int(a['TotalIntervalSize']) for a in amplicons.values()]),
                                sum([int(a['#BreakpointEdges']) for a in amplicons.values()]),
                                len([a for a in amplicons.values() if a['OncogenesAmplified'].strip() != ',']) != 
                                    0,
                                res[key_map['meta']]['clinical_data']['vital_status'],
                                relapse)
                                )
    else:
      output.write("%s,%s,%f,%d,%d,%d,%s,%s,%s\n" % (res[key_map['WGS__patient_id']],
                                disease_type,amp,  
                                res[key_map['meta']]['clinical_data']['days_to_last_known_alive'],
                                0,
                                0,
                                False,
                                res[key_map['meta']]['clinical_data']['vital_status'],
                                relapse)
                                )
  output.close()

def read_summary_file(summary_file, cnv_file=None):
  input = open(summary_file,'r')
  amplicons = {}
  line = input.next()
  line = input.next()
  id = None
  for line in input:
    res = line.strip().split(' ')
    if len(res) == 2:
      continue
    if res[1] == 'AmpliconID':
      id = int(res[3])
      amplicons.setdefault(id,{})      
    else:
      amplicons[id].setdefault(res[1],res[3])
  if len(amplicons[id].keys()) < 5:
    return None
  return amplicons

def check_amplicons():
  ensembl = load_ensemble_genes()
  files = os.listdir('%s/tmp/' % PANCANCER_DIR)
  files = [f for f in files if f.find("summary.txt") != -1]
  output = open('%s/analysis/amplicons.csv' % PANCANCER_DIR,'w')
  output.write('%s,%s,%s,%s,%s\n' % ('Sample','Amplicon','Interval','Type','Gene'))
  for f in files:
    sample = f[0:-12]
    lines = [line.strip().split(' ')[-1] for line in open('%s/tmp/%s' % (PANCANCER_DIR,f), 'r') if line.find('Intervals = ') != -1 & line.find('#') == -1]
    hg_ints = interval_list()
    idx = 0
    for line in lines:
      idx=idx+1
      for i in line.split(','):
        hg_int = hg19.interval(i)
        hg_int.chrom = "chr%s" % hg_int.chrom
        genes = Set([h[0].info['gene_id'] for h in ensembl.intersection([hg_int])])
        for g in genes:
          output.write('%s,%d,%s,%s,%s\n' % (sample,idx,"%s:%d-%d" % (hg_int.chrom, hg_int.start, hg_int.end),"Other",g))
  output.close()

def onco_ensemble():
  ensembl = load_ensemble_genes()
  oncos = {}
  for g in hg19.oncogene_list:
    hits = [oncos.setdefault(h[0].info['gene_id'],h[0]) for h in ensembl.intersection([g])]
  output = open('/pedigree2/projects/namphuon/programs/pancancer/analysis/onco_ensembl.csv','w')
  output.write('Gene\n')
  for g in oncos.keys():
    output.write(g+"\n")
  output.close()

def compute_fpkm_cnv():
  ensembl = load_ensemble_genes()
  meta_file = '/pedigree/projects/extrachromosome/data/TCGA/metadata/masked/masked_filemap.tsv'
  type = 'TCGA-GBM'  
  lines = [line.strip().split('\t') for line in open(meta_file)]
  lines = [line for line in lines if line[1] == 'Primary Tumor' and line[5] == type]
  meta = dict([(line[4][0:-4],line) for line in lines])
  samples = [s for s in os.listdir('/pedigree2/projects/namphuon/data/cancer_viral_integration/data/TCGA/TCGA-GBM/') if s in meta and os.path.exists('/pedigree2/projects/namphuon/data/cancer_viral_integration/data/TCGA/TCGA-GBM/%s/tumor/fpkm-uq/%s.txt' % (s,s))]
  data_map = {}
  #load up fpkm-uq values/cnv map
  for sample in samples: 
    print "Loading %s" % sample
    #load up fpkm-uq values  
    fpkm_file = '/pedigree2/projects/namphuon/data/cancer_viral_integration/data/TCGA/TCGA-GBM/%s/tumor/fpkm-uq/%s.txt' % (sample,sample)   
    data_map.setdefault(sample,{})['fpkm'] = load_fpkm(fpkm_file)
    
    #load up cnv_map  
    cnv_file = "/pedigree/projects/extrachromosome/data/TCGA/cnvsegs_hg38/masked/%s" % meta[sample][3].replace('.txt','.bed')
    data_map.setdefault(sample,{})['cnv'] = load_cnv_file(cnv_file)
    data_map[sample]['cnv'].sort()
  
  #Now for each sample and each gene, compute CNV versus FPKM
  ensembl_ids = dict([(g.info['gene_id'],g) for g in ensembl])
  genes = [g for g in data_map[sample ]['fpkm'].keys() if g in ensembl_ids]
  
  output = open('/pedigree2/projects/namphuon/programs/pancancer/analysis/tmp/cnv_fpkm.csv.%s' % sys.argv[1],'w')
  if sys.argv[1] == 0:
    output.write("%s,%s,%s,%s,%s,%s\n" % ('Sample','Type','Gene','FPKM_UQ','CNV','AA_Type'))
  idx = -1
  genes.sort()
  samples.sort()
  for g in genes:
    idx+=1
    if idx % 10 != int(sys.argv[1]):
      continue      
    for s in samples:
      hits = [c[0] for c in data_map[s]['cnv'].intersection([ensembl_ids[g]])]
      cnv = compute_average_amplification(ensembl_ids[g], hits)
      fpkm = data_map[s]['fpkm'][g]
      output.write("%s,%s,%s,%s,%s,%s\n" % (s,type,g,fpkm,cnv,'unknown'))
  output.close()

def compute_average_amplification(hg_int, cnv):
  avg = 0
  positions = 0
  for c in cnv:
    c_copy = hg19.interval(c.chrom,c.start,c.end,info=c.info)
    if c_copy.start < hg_int.start:
      c_copy.start = hg_int.start
    if c_copy.end > hg_int.end:
      c_copy.end = hg_int.end
    length = c_copy.end-c_copy.start
    avg+=length*(c_copy.info['Amplification'])
    positions+=length
  leftover = (hg_int.end-hg_int.start)-positions
  avg+=leftover
  return avg/(hg_int.end-hg_int.start)

def load_cnv_file(input):
  intervals = hg19.interval_list()
  for line in open(input):
    res = line.strip().split('\t')
    intervals.append(hg19.interval(res[0],res[1],res[2],info={'Amplification':2**float(res[4])}))
  intervals.sort()
  return intervals

def run_samples_python(file = '/home/namphuon/programs/pancancer/output/all_bam_cnv.csv', filter = 'TCGA-GBM', min_cnv = 3, min_size = 50000):
  lines = [line.strip().split('\t') for line in open(file)]
  counter = -1
  for line in lines[::-1]:
    if line[1] != filter:
      continue
    #Check to see if sample has been run
    output_dir = "output/TCGA/%s/%s/%s" % (line[1],line[0],'tumor')
    try:
      res = subprocess.check_output(['gsutil','ls', 'gs://aa-data-repo/%s/%s_summary.txt' % (output_dir,line[0])],stderr=subprocess.STDOUT)
      if len(res.strip()) != 0:
        continue
    except subprocess.CalledProcessError as error:
      print "Try to analyze %s" % line[0]
      #Analyze this sample 
    #Check to see if CNV file has been downloaded    
    counter+=1    
    if counter % 4 != int(sys.argv[1]):
      continue
    try:
      res = subprocess.check_output(['gsutil','cat', 'gs://aa-data-repo/TCGA/%s/%s/tumor/%s.hg19.%s_%s.masked_cnv.bed' % (line[1],line[0],line[0],min_cnv,min_size)],stderr=subprocess.STDOUT)
      if len(res.strip()) == 0:
        continue
    except subprocess.CalledProcessError as error:
      continue      
    path = tempfile.mkdtemp(dir=TEMP_DIR)  
    
    #Copy files from path
    print "Copying files" 
    try:
      res = subprocess.check_output(['gsutil','cp', 'gs://aa-data-repo/TCGA/%s/%s/tumor/%s.hg19.%s_%s.masked_cnv.bed' % (line[1],line[0],line[0],min_cnv,min_size), '%s/%s.bed' % (path,line[0])],stderr=subprocess.STDOUT)
      res = subprocess.check_output(['gsutil','cp', '%s.bai' % line[3], '%s/%s.bam.bai' % (path,line[0])],stderr=subprocess.STDOUT)
      res = subprocess.check_output(['gsutil','cp', '%s' % line[3], '%s/%s.bam' % (path,line[0])],stderr=subprocess.STDOUT)
      os.system('mkdir -p %s/%s/' % (path,line[0]))
      #os.system("python /home/namphuon/programs/docker_env/AmpliconArchitect/src/AmpliconArchitect.py --ref GRCh37 --downsample 10 --bam %s/%s.bam --bed %s/%s.bed --out %s/%s/%s > %s/%s/log.%s 2>&1" % (path,line[0],path,line[0],path,line[0],line[0], path, line[0], line[0]))
      os.system("python /pedigree2/projects/namphuon/programs/docker_env/AmpliconArchitect/src/AmpliconArchitect.py --ref GRCh37 --downsample 10 --bam %s/%s.bam --bed %s/%s.bed --out %s/%s/%s > %s/%s/log.%s 2>&1" % (path,line[0],path,line[0],path,line[0],line[0], path, line[0], line[0]))
      os.system('gsutil cp %s/%s/* gs://aa-data-repo/%s/' % (path,line[0],output_dir))
      os.system('rm -rf %s/' % path)
    except:
      print "Failed %s" " ".join(['gsutil','cp', '%s' % line[3], '%s/%s.bam' % (path,line[0])])
      os.system('rm -rf %s/' % path)    
      continue

def run_samples(file = '/home/namphuon/programs/pancancer/output/all_bam_cnv.csv', filter = 'TCGA-GBM', min_cnv = 3, min_size = 50000):
  lines = [line.strip().split('\t') for line in open(file)]
  path = tempfile.mkdtemp()
  output = open('%s/scripts/gbm.sh' % PANCANCER_DIR, 'w')
  for line in lines:
    if line[1] != filter:
      continue
    output_dir = "output/TCGA/%s/%s/%s" % (line[1],line[0],'tumor')      
    try:
      res = subprocess.check_output(['gsutil','ls', 'gs://aa-data-repo/%s/%s_summary.txt' % (output_dir,line[0])],stderr=subprocess.STDOUT)
      if len(res.strip()) != 0:
        continue
    except subprocess.CalledProcessError as error:
      print "Try to analyze %s" % line[0]
      #Analyze this sample 
    #Check to see if CNV file has been downloaded
    try:
      res = subprocess.check_output(['gsutil','cat', 'gs://aa-data-repo/TCGA/%s/%s/tumor/%s.hg19.%s_%s.masked_cnv.bed' % (line[1],line[0],line[0],min_cnv,min_size)],stderr=subprocess.STDOUT)
      if len(res.strip()) == 0:
        continue
    except subprocess.CalledProcessError as error:
      continue       
    res = subprocess.check_output(['gsutil','ls', '-lh', line[3]],stderr=subprocess.STDOUT).strip()
    mem = int(float(res.split(' ')[0])+25)
    output.write("dsub --image us.gcr.io/aa-test-175718/aa --preemptible --project aa-test-175718 --zones \"us-west1-*\" --logging gs://aa-data-repo/%s/logging/ \\\n\t--output LOG_FILE=gs://aa-data-repo/%s/log \\\n\t --input-recursive DATA_REPO=gs://aa-data-repo/data_repo/ \\\n\t--output-recursive OUTPUT_DIR=gs://aa-data-repo/%s --disk-size %d \\\n\t--input SH=%s BAM_FILE=%s BAI_FILE=%s.bai BED_FILE=gs://aa-data-repo/TCGA/%s/%s/tumor/%s.hg19.%s_%s.masked_cnv.bed \\\n\t--env PREFIX=%s --command 'sh ${SH}'\n" %  (output_dir,output_dir,output_dir, mem,"gs://aa-data-repo/data_repo/test.sh",line[3],line[3],line[1],line[0],line[0],min_cnv,min_size,line[0]))
  output.close()  

def prepare_samples(file = '/home/namphuon/programs/pancancer/output/all_bam_cnv.csv', filter = 'TCGA-GBM', min_cnv = 3, min_size = 50000):
  lines = [line.strip().split('\t') for line in open(file)]
  for line in lines:
    path = tempfile.mkdtemp()  
    if line[1] != filter:
      continue
    #Check to see if CNV file has been downloaded
    try:
      res = subprocess.check_output(['gsutil','ls', 'gs://aa-data-repo/TCGA/%s/%s/tumor/%s.masked_cnv.txt' % (line[1],line[0],line[0])],stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:      
      #Download the file
      os.system('gdc-client download -d %s %s' % (path, line[-1]))
      t = "%s/%s" % (path, os.listdir(path)[0])
      cnv_file = [f for f in os.listdir(t) if f.find('seg.txt') != -1][0]
      subprocess.check_output(['gsutil','cp', '%s/%s' % (t,cnv_file), 'gs://aa-data-repo/TCGA/%s/%s/tumor/%s.masked_cnv.txt' % (line[1],line[0],line[0])],stderr=subprocess.STDOUT)    
    try:
      res = subprocess.check_output(['gsutil','ls', 'gs://aa-data-repo/TCGA/%s/%s/tumor/%s.masked_cnv.txt' % (line[1],line[0],line[0])],stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:      
      "Failed to download gs://aa-data-repo/TCGA/%s/%s/tumor/%s.masked_cnv.txt" % (line[1],line[0],line[0]) 
      continue    
    try:
      res = subprocess.check_output(['gsutil','ls', 'gs://aa-data-repo/TCGA/%s/%s/tumor/%s.hg19.%s_%s.masked_cnv.bed' % (line[1],line[0],line[0],min_cnv,min_size)],stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
      tfile = "%s/%s.tmp" % (path,line[0])
      res = subprocess.check_output(['gsutil','cp', 'gs://aa-data-repo/TCGA/%s/%s/tumor/%s.masked_cnv.txt' % (line[1],line[0],line[0]), tfile],stderr=subprocess.STDOUT)
      os.system("awk '{print $2,$3,$4,$5,$6}' %s | tail -n +2 > %s.bed" % (tfile, tfile))
      ofile = '%s/%s.hg19.%s_%s.masked_cnv.bed' % (path,line[0],min_cnv,min_size)
      find_amplified_intervals("%s.bed" % (tfile), ofile, min_cnv=min_cnv, min_size=min_size)
      subprocess.check_output(['gsutil','cp', '%s' % (ofile), 'gs://aa-data-repo/TCGA/%s/%s/tumor/' % (line[1],line[0])],stderr=subprocess.STDOUT) 
    os.system('rm -rf %s' % path)

def find_amplified_intervals(input,output,min_cnv=5,min_size=10000):
  out = open(output,'w')
  for line in open(input):
    res = line.strip().split(' ')
    if 2**float(res[4]) >= min_cnv and (int(res[2])-int(res[1]) >= min_size):
      out.write(line)
  out.close()

def collect_all_rnaseq_files():
  files_endpt = 'https://api.gdc.cancer.gov/files'
  filt = {"op":"and",
            "content":[
            {"op":"=",
              "content":{
                  "field": "data_type",
                  "value": "Gene Expression Quantification"
              }}
            ,
            {"op":"or",
              "content":[
              {"op":"=",
                "content":{
                    "field": "cases.samples.sample_type",
                    "value": "*umor*"
                }},
              {"op":"=",
                "content":{
                    "field": "cases.samples.sample_type",
                    "value": "*ancer*"
                }},                
              {"op":"=",
                "content":{
                    "field": "cases.samples.sample_type",
                    "value": "*static*"
                }}]              
            }
            ]}                                 
  params = {'filters':json.dumps(filt),'size':100000, 'fields':"cases.project.project_id,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,file_id,file_name,cases.submitter_id,cases.samples.portions.portion_id,cases.samples.annotations.annotation_id,cases.samples.annotations.case_id,cases.samples.annotations.entity_id,cases.samples.portions.submitter_id,analysis.workflow_type,analysis.submitter_id,analysis.analysis_id,analysis.input_files.submitter_id,cases.submitter_aliquot_ids,cases.annotations.entity_id"}
  response = requests.get(files_endpt, params = params)  
  return response

def merge_picked_results_rnaseq_fpkm():
  sample_mapper = {}
  rs = [r for r in response.json()['data']['hits']]
  for r in response.json()['data']['hits']:
    if r['analysis']['workflow_type'] != 'HTSeq - FPKM-UQ' or r['cases'][0]['project']['project_id'].find('TCGA') == -1:
      continue    
    portion_id = r['cases'][0]['samples'][0]['portions'][0]['submitter_id']
    if portion_id in sample_mapper:
      print "Fuck the TCGA"
    sample_mapper[portion_id] = r
  for c in complete.values():
    if c[key_map['WGS__portion_id']] in sample_mapper:
      c.append(sample_mapper[c[key_map['WGS__portion_id']]]['file_id'])
  keys_to_remove = []
  for (k,v) in complete.items():
    if len(v) == 87:
      print c[key_map['rnaseq_fpkm_uq']]
    else:
      keys_to_remove.append(k)
      print "Missing"
  for c in keys_to_remove:
    del complete[c]

def collect_all_metadata_files():
  files_endpt = 'https://api.gdc.cancer.gov/files'
  filt = {"op":"and",
            "content":[
            {"op":"=",
              "content":{
                  "field": "data_type",
                  "value": "Masked Copy Number Segment"
              }}
            ,
            {"op":"or",
              "content":[
              {"op":"=",
                "content":{
                    "field": "cases.samples.sample_type",
                    "value": "*umor*"
                }},
              {"op":"=",
                "content":{
                    "field": "cases.samples.sample_type",
                    "value": "*ancer*"
                }},                
              {"op":"=",
                "content":{
                    "field": "cases.samples.sample_type",
                    "value": "*static*"
                }}]              
            }
            ]}                                 
  params = {'filters':json.dumps(filt),'size':100000, 'fields':"cases.project.project_id,cases.samples.sample_type,analysis.workflow_type,cases.samples.submitter_id,cases.samples.sample_id,file_id,file_name,cases.submitter_id"}
  # requests URL-encodes automatically
  
  response = requests.get(files_endpt, params = params)  
  output = open('/home/namphuon/data/dbGap/metadata/masked.cnv.csv','w')
  for hit in response.json()['data']['hits']:
    output.write('%s\t%s\t%s\t%s\t%s\n' % (hit['cases'][0]['project']['project_id'], hit['cases'][0]['submitter_id'], hit['cases'][0]['samples'][0]['sample_type'], hit['id'],hit['file_id']))
  output.close()  
  
  filt = {"op":"and",
            "content":[
            {"op":"=",
              "content":{
                  "field": "data_type",
                  "value": "Aligned reads"
              }}
            ,
            {"op":"=",
              "content":{
                  "field": "data_format",
                  "value": "BAM"
              }}
            ,
            {"op":"or",
              "content":[
                {"op":"=",
                "content":{
                    "field": "experimental_strategy",
                    "value": "WGS"
                }},
                {"op":"=",
                "content":{
                    "field": "experimental_strategy",
                    "value": "RNA-Seq"
                }}]
            }
            ,
            {"op":"or",
              "content":[
              {"op":"=",
                "content":{
                    "field": "cases.samples.sample_type",
                    "value": "*umor*"
                }},
              {"op":"=",
                "content":{
                    "field": "cases.samples.sample_type",
                    "value": "*ancer*"
                }},                
              {"op":"=",
                "content":{
                    "field": "cases.samples.sample_type",
                    "value": "*static*"
                }}]              
            }
            ]}                                 
  params = {'filters':json.dumps(filt),'size':1000000, 'fields':"cases.project.project_id,cases.samples.sample_type,cases.samples.submitter_id,cases.samples.sample_id,file_id,file_name,cases.submitter_id,experimental_strategy,data_format,data_type"}
  # requests URL-encodes automatically
  files_endpt = 'https://api.gdc.cancer.gov/legacy/files'  
  response = requests.get(files_endpt, params = params)  
  output = open('/home/namphuon/data/dbGap/metadata/bam.csv','w')
  for hit in response.json()['data']['hits']:
    output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (hit['cases'][0]['project']['project_id'], hit['cases'][0]['submitter_id'], hit['cases'][0]['samples'][0]['sample_type'], hit['id'],hit['file_id'],hit['experimental_strategy'],hit['file_name']))
  output.close()

def get_legacy_samples(file='/home/namphuon/data/dbGap/legacy.bam.json'):
  bams = json.load(open(bam_json_file))  
  samples = Set()
  for bam in bams:      
    samples.add(bam['cases'][0]['submitter_id'])
  return samples

def prepare_tcga_data(bam_csv="/home/namphuon/data/dbGap/metadata/bam.csv", cnv_csv="/home/namphuon/data/dbGap/metadata/masked.cnv.csv"):
  samples = {}
  lines =[line.strip().split('\t') for line in open(bam_csv,'r')]  
  for line in lines:
    samples.setdefault(line[1],{}).setdefault(line[2],{}).setdefault(line[5],line)
  
  lines =[line.strip().split('\t') for line in open(cnv_csv,'r')]  
  for line in lines:
    samples.setdefault(line[1],{}).setdefault(line[2],{}).setdefault("Masked CNV",line)
  
  for sample in samples.keys():    
    for type in samples[sample].keys():
      if type != 'Primary Tumor':
        del samples[sample][type]
        continue
      if 'WGS' not in samples[sample][type] or 'RNA-Seq' not in samples[sample][type] or 'Masked CNV' not in samples[sample][type]:
        del samples[sample][type]
    if len(samples[sample].keys()) == 0:
      del samples[sample]        
  type = 'Primary Tumor'  
  strategies = ['WGS', 'RNA-Seq']  
  for sample in samples:
    for t in strategies:
      if samples[sample][type][t][-1].find('gs') != -1:    
        continue
      for key in ['A','B']:
        data = service.samples().cloud_storage_file_paths(sample_barcode='%s-01%s' % (sample,key), experimental_strategy=t).execute()      
        if data['count'] == 0:
          continue
        gs = [g for g in data['cloud_storage_file_paths'] if g.find(samples[sample][type][t][-1]) != -1]
        if len(gs) == 0:
          continue
        else:
          samples[sample][type][t].append(gs[0])        
          break
      if samples[sample][type][t][-1].find('gs') == -1:
        print 'Failed to find gs for %s %s' % (sample,t)
  for sample in samples.keys():
    for t in strategies:
      if samples[sample][type][t][-1].find('gs') == -1:    
        del samples[sample]
        break
  output = open('%s/output/all_bam_cnv.csv' % PANCANCER_DIR, 'w')
  type = 'Primary Tumor'  
  strategies = ['WGS', 'RNA-Seq','Genotyping Array']
  for s in samples.keys():    
    output.write("%s\t%s\t%s\t%s\n" % (s,samples[s][type]['WGS'][0],type,"\t".join([samples[s][type]['Masked CNV'][-1] if t == 'Genotyping Array' else samples[s][type][t][-1]  for t in strategies])))
  output.close()

def load_meta_data(bam_json_file, cnv_json_file):
  bams = json.load(open(bam_json_file))
  samples = {}
  for bam in bams:
    project = bam['cases'][0]['project']
    sample_type = bam['cases'][0]['samples'][0]['sample_type']
    sample_type_id = bam['cases'][0]['samples'][0]['sample_type_id']
    strategy = bam['experimental_strategy']
    file_id = bam['file_id']
    file_name = bam['file_name']
    sample_name = bam['cases'][0]['submitter_id']
    samples.setdefault(sample_name,{}).setdefault(sample_type,{}).setdefault(strategy, {'file_id':file_id,'project':project,'file_name':file_name})
  
  cnvs = json.load(open(cnv_json_file))
  for cnv in cnvs:
    project = cnv['cases'][0]['project']
    sample_type = cnv['cases'][0]['samples'][0]['sample_type']
    sample_type_id = cnv['cases'][0]['samples'][0]['sample_type_id']
    strategy = cnv['experimental_strategy']
    file_id = cnv['file_id']
    file_name = cnv['file_name']
    sample_name = cnv['cases'][0]['submitter_id']
    samples.setdefault(sample_name,{}).setdefault(sample_type,{}).setdefault(strategy, {'file_id':file_id,'project':project,'file_name':file_name})
  return samples

def filter_samples(samples):
  type = 'Primary Tumor'
  strategies = ['Genotyping Array', 'RNA-Seq', 'WGS']
  keep = [s for s in samples if len([st for st in samples[s][type] if st in strategies]) == len(strategies)]
  for s in samples.keys():
    if s not in keep:
      del samples[s]
  return samples

def write_metadata(samples,file='%s/output/hnsc.csv' % PANCANCER_DIR):
  output = open(file, 'w')
  type = 'Primary Tumor'  
  strategies = ['WGS', 'RNA-Seq','Genotyping Array']
  for s in samples.keys():
    output.write("%s\t%s\t%s\n" % (s,type,"\t".join([samples[s][type][t]['file_id'] if t == 'Genotyping Array' else samples[s][type][t]['gs']  for t in strategies])))
  output.close()

def pull_metadata_files(samples):
  type = 'Primary Tumor'  
  strategies = ['WGS', 'RNA-Seq']
  
  for sample in samples:
    for t in strategies:
      for key in ['A','B']:
        data = service.samples().cloud_storage_file_paths(sample_barcode='%s-01%s' % (sample,key), experimental_strategy=t).execute()      
        if data['count'] == 0:
          continue
        gs = [g for g in data['cloud_storage_file_paths'] if g.find(samples[sample][type][t]['file_name']) != -1]
        if len(gs) != 1:
          continue
        else:
          samples[sample][type][t]['gs'] = gs[0]        
          break
      if 'gs' not in samples[sample][type][t]:
        print 'Failed to find gs for %s %s' % (sample,t)

def get_unauthorized_service(api = 'isb_cgc_tcga_api'):
    version = 'v3'
    site = "https://api-dot-isb-cgc.appspot.com"
    discovery_url = '%s/_ah/api/discovery/v1/apis/%s/%s/rest' % (site, api, version)
    return build(api, version, discoveryServiceUrl=discovery_url, http=httplib2.Http())

def load_ensemble_genes():
  #file = '/pedigree2/projects/namphuon/data/references/hg19/annotations/Homo_sapiens.GRCh38.89.chr.gtf'
  file = '/pedigree2/projects/namphuon/data/references/hg19/annotations/Homo_sapiens.GRCh38.89.chr.gtf'
  lines = [line.strip().split('\t') for line in open(file,'r') if len(line.strip().split('\t')) > 4 and line.strip().split('\t')[2] == 'gene']
  ensemble_gene_list = hg19.interval_list()  
  foo = [ensemble_gene_list.append(hg19.interval("chr%s" % line[0], int(line[3]), int(line[4]),strand=1 if line[6]=='+' else -1, info={'data':line})) for line in lines]
  foo = [g.info.setdefault(s.strip().split(' ')[0],s.strip().split(' ')[1].replace('"',"")) for g in ensemble_gene_list for s in g.info['data'][8][:-1].split(';') ]  
  ensemble_gene_list.sort()
  return ensemble_gene_list

def load_fpkm(file,names=Set()):
  input = open(file, 'r')
  lines = [line.strip().split('\t') for line in input]
  if len(names) != 0:
    fpkm = dict([(res[0].split('.')[0],float(res[1])) for res in lines if res[0].split('.')[0] in names])
  else:
    fpkm = dict([(res[0].split('.')[0],float(res[1])) for res in lines])
  return fpkm

def prefix_chr(opened_bam_file):
  chr_count = [i['SN'] for i in opened_bam_file.header['SQ'] if 'SN' in i and i['SN'].find('chr') != -1]
  return len(chr_count) > 20  

def read_read_counts(index):
  read_counts = dict()
  ks = complete.keys()
  #ks = [k for k in complete.keys() if complete[k][key_map['WGS__cases__project__project_id']] == disease]
  for k in ks:
    res=complete[k]
    tumor_type = res[key_map['WGS__sampleType']]
    disease_type = res[key_map['WGS__cases__project__project_id']]
    out_dir = '%s/analyses/tumor/%s/min_cnv1/%s/'  % (PANCANCER_DIR, disease_type, res[key_map['outdir_basename']])          
    if tumor_type != 'TP':
      continue
    flagstat_file = "%s/flagstat.csv" % (out_dir)
    read_file = "%s/read_counts.csv" % (out_dir)
    fpkm_file = "%s/fpkm_uq.csv" % (out_dir)
    if not os.path.exists(flagstat_file) or not os.path.exists(read_file) or not os.path.exists(fpkm_file):
      continue
    flagstat = read_samtools_flagstat(flagstat_file)
    reads = read_readcounts(read_file)
    read_counts[k] = (reads, flagstat)            
    cnv_bed_file =  '%s/%s.hg19.bed' % (out_dir,res[key_map['WGS__patient_id']])
    cnv = load_cnv_file(cnv_bed_file)
    fpkm = load_fpkm(fpkm_file)
    all_map = load_amplicon_file(res)
    read_counts[k] = (reads, flagstat, cnv, fpkm, all_map)
  #Compute FPKM for any segment of a cycle that is cyclic
  output = open('%s/analyses/gene_amplicon_tpm.%d.csv' % (PANCANCER_DIR, index),'w')
  #output.write('Sample,Disease,Gene_name ,Gene,Structure,FPKM,CNV,Oncogene\n')  
  counter = -1
  for disease in diseases:
    counter+=1
    if counter % 6 != index:
      continue
    cycle_maps = {}
    segment_maps = {}
    interval_maps = {}
    samples = [v[key_map['Tag']] for v in complete.values() if v[key_map['WGS__cases__project__project_id']] == disease and v[key_map['Tag']] in read_counts]    
    for reference_sample in samples:
      for (reference_amp_id,all_map) in read_counts[reference_sample][4].items():
        cycles = [(cid,cmap) for (cid,cmap) in all_map['cycle_map'].items() if cmap['cycle'][0][0] != '0' and cmap['cycle'][-1][0] != '0']
        interval_maps.setdefault(reference_sample, hg19.interval_list()).extend(all_map['interval_map'].values())
        interval_maps[reference_sample].sort()
        if len(cycles) == 0:
          continue
        for c in cycles:
          cycle_interval_list = hg19.interval_list()        
          for segment in c[1]['cycle']:
            if segment[:-1] == '0':
              continue
            c_seg = all_map['segment_map'][segment[:-1]]
            c_seg.info['cycle_id'] = c[0]
            c_seg.info['cycle_map'] = c[1]            
            cycle_interval_list.append(c_seg)
          cycle_interval_list.sort()
#         if sum([(c.end-c.start) for c in cycle_interval_list]) > 50000:
#           continue
          cycle_maps.setdefault(reference_sample,{}).setdefault(reference_amp_id,{}).setdefault(
                                c[0], hg19.interval_list()).extend(cycle_interval_list)
          segment_maps.setdefault(reference_sample,hg19.interval_list()).extend(cycle_interval_list)
          segment_maps[reference_sample].sort()
    #Find all intersecting genes
    all_intervals = hg19.interval_list([i for v in interval_maps.values() for i in v])
    all_intervals.sort()
    gene_list = [h[0] for h in ensembl.intersection(all_intervals)]
    gene_map = {}
    foo = [gene_map.setdefault(g.info['gene_id'],hg19.interval_list()).append(g) for g in gene_list]
    for (k,v) in gene_map.items():
      gene_map[k] = hg19.interval_list(Set(gene_map[k]))
      gene_map[k].sort()
      if len(gene_map[k]) > 1:
        break
    for (gene_id, g) in gene_map.items():
      oncogene = len(hg19.oncogene_list.intersection(g)) > 0
      for current_sample in samples:        
        if g[0].info['gene_id'] not in read_counts[current_sample][3]:        
          continue      
        structure_type = "None"
        if current_sample in cycle_maps and len(segment_maps[current_sample].intersection(g)) > 0:
          structure_type = "Cycle"
        elif current_sample in interval_maps and len(interval_maps[current_sample].intersection(g)) > 0:
          structure_type = "Amplicon"
        cnv_hits = hg19.interval_list(Set([c[0] for c in hg19.interval_list(read_counts[current_sample][2].intersection(g))]))
        cnv_hits.sort()
        
        cnv = compute_average_amplification(g[0], cnv_hits)
        fpkm = read_counts[current_sample][3][g[0].info['gene_id']]
        output.write('%s,%s,%s,%s,%s,%f,%f,%s\n' % (current_sample, disease, 
        g[0].info['gene_name'], g[0].info['gene_id'], structure_type, fpkm, cnv, oncogene))
  output.close()    


service = get_unauthorized_service(api='isb_cgc_tcga_api')
(data_map, key_map) = load_meta_file()
#(complete, keeper) = get_rnaseq_bucket(data_map, service, key_map)
#pickle.dump( complete, open( "%s/analyses/complete.p" % PANCANCER_DIR, "wb" ) )
#pickle.dump( keeper, open( "%s/analyses/keeper.p" % PANCANCER_DIR, "wb" ) )
complete = pickle.load( open(  "%s/analyses/complete.p" % PANCANCER_DIR, "rb" ) )
keeper = pickle.load( open(  "%s/analyses/keeper.p" % PANCANCER_DIR, "rb" ) )
#download_fpkm_values(complete,key_map)
ensembl = load_ensemble_genes()
for h in hg19.oncogene_list:
  h.chrom = "chr%s" % h.chrom
read_read_counts(index=int(sys.argv[1]))
