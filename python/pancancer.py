from __future__ import division
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

def load_meta_file():
  data_map = dict()
  #for disease in [d for d in os.listdir('%s/analyses/tumor' % PANCANCER_DIR) if d.find('TCGA') != -1]:
  for disease in ['TCGA-GBM', 'TCGA-LGG']:
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
    for i in input:
      res = [r.replace('"','') for r in i.strip().split('\t')]
      data_map.setdefault(res[0],res)
    #Find if the file has an RNA-seq file
  return (data_map, key_map)
  
def build_amplicon_list(complete):
  for disease in [d for d in os.listdir('%s/analyses/tumor' % PANCANCER_DIR) if d.find('TCGA') != -1]:
    collect_amplicons(complete,disease=disease)
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
    output.write("dsub --image us.gcr.io/aa-test-175718/read_count --preemptible --project aa-test-175718 --zones \"us-west1-*\" --logging gs://aa-data-repo/%s/logging/ \\\n\t--output LOG_FILE=gs://aa-data-repo/%s/log.txt \\\n\t \\\n\t--output-recursive OUTPUT_DIR=gs://aa-data-repo/%s --disk-size %d \\\n\t--input SEGMENT_FILE=gs://aa-data-repo/metadata/amplicons.%s.csv BAM_FILE=%s BAI_FILE=%s.bai \\\n\t --command 'python /home/read_count.py -b ${BAM_FILE} -s ${SEGMENT_FILE} -f ${OUTPUT_DIR}/flagstat.csv -o ${OUTPUT_DIR}/read_counts.csv > ${LOG_FILE} 2>&1; samtools flagstat ${BAM_FILE} > ${OUTPUT_DIR}/flagstat.csv'\n" %  (out_dir,out_dir,out_dir, mem, disease, rna_bam_file,rna_bam_file))
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
      keeper[d[0]] = tmp
    else:
      tmp = copy.copy(d)
      tmp.append(None)
      tmp.append(res)      
      tmp.append(cases)            
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

def collect_amplicons(complete,disease='TCGA-GBM'):
  amplicons = hg19.interval_list()
  for k in complete.keys():
    res=complete[k]
    tumor_type = res[key_map['WGS__sampleType']]
    disease_type = res[key_map['WGS__cases__project__project_id']]
    out_dir = '%s/analyses/tumor/%s/min_cnv1/%s/'  % (PANCANCER_DIR, disease_type, res[key_map['outdir_basename']])          
    if tumor_type != 'TP' or disease_type != disease:
      continue
    summary_file = "%s/output/%s-MINCNV1_summary.txt" % (out_dir,res[key_map['Tag']])    
    if os.path.exists(summary_file):    
      amps = read_summary_file(summary_file)
    for (ai, a) in amps.items():
      ints = a['Intervals'].split(',')
      for i in ints:
        amplicons.append(hg19.interval(i,info={'amplicon':a, 'amplicon_id':ai, 'meta':res}))
  amplicons.sort()
  amplicons = hg19.interval_list([a for a in amplicons if a.chrom in chrs])
  
  output = open('%s/metadata/amplicons.%s.csv' % (PANCANCER_DIR,disease),'w')
  output.write('Sample,Disease,Amplicon_id,Segment\n')
  for a in amplicons:
    id = a.info['meta'][0]
    disease = a.info['meta'][key_map['WGS__cases__project__project_id']]
    amp = a.info['amplicon_id']
    seg = "chr%s:%d-%d" % (a.chrom, a.start, a.end)
    output.write('%s,%s,%s,%s\n' % (id, disease, amp, seg))
  output.close()
  
  output = open('%s/metadata/complete.%s.csv' % (PANCANCER_DIR,disease),'w')
  output.write('Sample,Disease,Bucket\n')
  for k in complete.keys():
    res=complete[k]
    id = res[0]
    tumor_type = res[key_map['WGS__sampleType']]
    disease_type = res[key_map['WGS__cases__project__project_id']]
    rna_bam_file = res[key_map['rnaseq_gsc']]
    if tumor_type != 'TP':
      continue    
    output.write('%s,%s,%s\n' % (id, disease, rna_bam_file))
  output.close()

def read_read_counts(complete, amplicons, disease='TCGA-GBM'):
  read_counts = dict()
  for k in complete.keys():
    res=complete[k]
    tumor_type = res[key_map['WGS__sampleType']]
    disease_type = res[key_map['WGS__cases__project__project_id']]
    out_dir = '%s/analyses/tumor/%s/min_cnv1/%s/'  % (PANCANCER_DIR, disease_type, res[key_map['outdir_basename']])          
    if tumor_type != 'TP':
      continue
    flagstat_file = "%s/flagstat.csv" % (out_dir)
    read_file = "%s/read_counts.csv" % (out_dir)
    if not os.path.exists(flagstat_file) or not os.path.exists(read_file):
      continue
    flagstat = read_samtools_flagstat(flagstat_file)
    reads = read_readcounts(read_file)    
    read_counts[k] = (reads, flagstat)            
    cnv_bed_file =  '%s/%s.hg19.bed' % (out_dir,res[key_map['WGS__patient_id']])
    cnv = load_cnv_file(cnv_bed_file)
    read_counts[k] = (reads, flagstat, cnv)    
    
  #Now for each amplicon, determine which samples intersect it so baselines can be determined
  amps = [a for a in amplicons if a.info['meta'][key_map['WGS__cases__project__project_id']] == disease]
  comps = dict([(i,v) for (i,v) in complete.items() if v[key_map['WGS__cases__project__project_id']] == disease])
  sample_amps = {}
  [sample_amps.setdefault(a.info['meta'][0],hg19.interval_list()).append(a) for a in amps]
  [v.sort() for v in sample_amps.values()]
  #Split into samples with and without interval in amplicon for each amplicon, find baseline, compute fold change
  for a in amps:
    samples_not_amplicon = [sample for (sample,sample_list) in sample_amps.items() if len(sample_list.intersection([a])) == 0]
      
    #Av. read from amplicon region in non-amp samples
    idx = 0
#    for s in samples_not_amplicon:
#      total_mapped_read = read_counts[s][1]['flag

def read_readcounts(read_file):
  read_count_intervals = {}
  input = open(read_file, 'r')
  for line in input:
    res = line.strip().split(',')  
    foo = read_count_intervals.setdefault("%s_%s" % (res[0],res[2]),{}).setdefault(res[1],res)
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
      output.write("%s,%s,%f,%d,%d,%d,%s,%s,%s\n" % (res[key_map['WGS__patient_id']],
                                disease_type,amp,  
                                res[key_map['meta']]['clinical_data']['days_to_last_known_alive'],
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

def read_summary_file(summary_file):
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
  avg+=2*leftover
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

#service = get_unauthorized_service(api='isb_cgc_tcga_api')
#(data_map, key_map) = load_meta_file()
#(complete, keeper) = get_rnaseq_bucket(data_map, service, key_map)
#run_read_count(complete, disease='TCGA-GBM')
#run_samples_python()
#run_samples_python(file='/pedigree2/projects/namphuon/bin/perl/scripts/all_bam_cnv.csv')
#compute_fpkm_cnv()

#   rna_bam_file = "gs://5aa919de-0aa0-43ec-9ec3-288481102b6d/tcga/GBM/RNA/RNA-Seq/UNC-LCCC/ILLUMINA/UNCID_1534748.d8164b02-4b3c-454d-945b-2838edb1b5b1.sorted_genome_alignments.bam"
#   rna_bam = pysam.Samfile(rna_bam_file, 'rb')      
