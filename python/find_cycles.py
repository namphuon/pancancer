import hg19util as hg19
import sys, os

#Reads a summary file
def read_summary_file(summary_file):
  input = open(summary_file,'r')
  amplicons = {}
  line = next(input, None)
  if line is None:
    return None
  line = next(input, None)
  if line is None:
    return None
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

#Reads a cycle file
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
      #Sometimes a source/sink cycle will not start/end with 0, but have the 0 within the cycle
      #Rotate the result such that 0s always start/end the cycle in those cases
      if cycle[0][0] != '0' and len([c for c in cycle if c[0] == '0']) > 0:
        idx = cycle.index([c for c in cycle if c == '0+'][0])
        cycle = cycle[idx:] + cycle[:idx]
      cycle_map[cycle_name] = {'cycle':cycle,'copy_count':copy_count}
  input.close()
  return (segment_map, interval_map, cycle_map)    

#Loads the summary file and amplicons for an AA run
#prefix is the prefix of the summary/amplicon files and in_dir is the
#directory of the AA results
def load_aa_result(prefix, in_dir):
  summary_file = "%s/%s_summary.txt" % (in_dir, prefix)
  all_map = {}
  amps = read_summary_file(summary_file)
  
  for (ai, a) in amps.items():
    ints = a['Intervals'].split(',')
    cycle_file = "%s/%s_amplicon%s_cycles.txt" % (in_dir,prefix,ai)
    (segment_map, interval_map, cycle_map) = read_cycle_file(cycle_file)
    all_map.setdefault('amplicons',{}).setdefault(ai,{})['segment_map'] = segment_map
    all_map.setdefault('amplicons',{}).setdefault(ai,{})['interval_map'] = interval_map
    all_map.setdefault('amplicons',{}).setdefault(ai,{})['cycle_map'] = cycle_map
  all_map['summary'] = amps
  return all_map

#Returns the cyclic paths of an amplicon of at least 10kb
def get_cyclic_path(amplicon, threshold = 10000):
  paths = hg19.interval_list()
  for cycle in amplicon['cycle_map'].keys():
    length = sum([amplicon['segment_map'][s[0:-1]].end-amplicon['segment_map'][s[0:-1]].start for s in amplicon['cycle_map'][cycle]['cycle'] if s[0] != '0'])
    if (amplicon['cycle_map'][cycle]['cycle'][0][0] != '0' and length >= threshold):
      paths.extend([amplicon['segment_map'][s[0]] for s in amplicon['cycle_map'][cycle]['cycle']])
  return paths

#Returns whether an amplicon has a cycle of at least 10kb
def valid_cycle(amplicon, threshold = 10000):
  is_cyclic = False
  for cycle in amplicon['cycle_map'].keys():
    length = sum([amplicon['segment_map'][s[0:-1]].end-amplicon['segment_map'][s[0:-1]].start for s in amplicon['cycle_map'][cycle]['cycle'] if s[0] != '0'])
    is_cyclic = is_cyclic or (amplicon['cycle_map'][cycle]['cycle'][0][0] != '0' and length >= threshold)      
  return is_cyclic

aa_result = load_aa_result(in_dir='/pedigree2/projects/namphuon/programs/pancancer/analyses/final/TUMOR-WGS-GWSNP6NOCNV/RAW/TCGA-06-0125-01A-01_WGS335134030397_CNV19272_mincnv1_20170905-133947/output', prefix='TCGA-06-0125-01A-01_WGS335134030397_CNV19272-MINCNV1')

#Now find which amplicons have a cyclic cycle
for amplicon in aa_result['amplicons'].keys():
  print "%d %s" % (amplicon, valid_cycle(aa_result['amplicons'][amplicon]))