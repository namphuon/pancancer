from __future__ import division
import sys, os, subprocess
import hg19util as hg19
import json,requests
import argparse
import tempfile,shutil
import glob
import os,random,tempfile,shutil
import os.path
import copy,math,bisect
from sets import Set
from collections import Counter
from argparse import ArgumentParser
from googleapiclient.discovery import build
import httplib2
import pprint


PANCANCER_DIR=os.environ['PANCANCER']
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

def load_cnv_file(input):
  intervals = hg19.interval_list()
  for line in open(input):
    res = line.strip().split('\t')
    intervals.append(hg19.interval(res[0],res[1],res[2],info={'Amplification':2**float(res[4])}))
  return intervals

def run_samples_python(file = '/home/namphuon/programs/pancancer/output/all_bam_cnv.csv', filter = 'TCGA-GBM', min_cnv = 3, min_size = 50000):
  lines = [line.strip().split('\t') for line in open(file)]
  for line in lines:
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
    try:
      res = subprocess.check_output(['gsutil','cat', 'gs://aa-data-repo/TCGA/%s/%s/tumor/%s.hg19.%s_%s.masked_cnv.bed' % (line[1],line[0],line[0],min_cnv,min_size)],stderr=subprocess.STDOUT)
      if len(res.strip()) == 0:
        continue
    except subprocess.CalledProcessError as error:
      continue         
    path = tempfile.mkdtemp()  
    
    #Copy files from path
    print "Copying files"    
    try:
      res = subprocess.check_output(['gsutil','cp', 'gs://aa-data-repo/TCGA/%s/%s/tumor/%s.hg19.%s_%s.masked_cnv.bed' % (line[1],line[0],line[0],min_cnv,min_size), '%s/%s.bed' % (path,line[0])],stderr=subprocess.STDOUT)
      res = subprocess.check_output(['gsutil','cp', '%s.bai' % line[3], '%s/%s.bam.bai' % (path,line[0])],stderr=subprocess.STDOUT)
      res = subprocess.check_output(['gsutil','cp', '%s' % line[3], '%s/%s.bam' % (path,line[0])],stderr=subprocess.STDOUT)
      os.system('mkdir -p %s/%s/' % (path,line[0]))
      os.system("python /home/namphuon/programs/docker_env/AmpliconArchitect/src/AmpliconArchitect.py --downsample 10 --bam %s/%s.bam --bed %s/%s.bed --out %s/%s/%s > %s/%s/log.%s 2>&1" % (path,line[0],path,line[0],path,line[0],line[0], path, line[0], line[0]))
      os.system('gsutil cp %s/%s/* gs://aa-data-repo/%s/' % (path,line[0],output_dir))
      os.system('rm -rf %s/' % path)
    except:
      os.system('rm -rf %s/' % path)    
      continue
  
def run_samples(file = '/home/namphuon/programs/pancancer/output/all_bam_cnv.csv', filter = 'TCGA-GBM', min_cnv = 3, min_size = 50000):
  lines = [line.strip().split('\t') for line in open(file)]
  path = tempfile.mkdtemp()
  output = open('%s/scripts/gbm.sh' % PANCANCER_DIR, 'w')
  for line in lines:
    if line[1] != filter:
      continue
    #Check to see if CNV file has been downloaded
    try:
      res = subprocess.check_output(['gsutil','cat', 'gs://aa-data-repo/TCGA/%s/%s/tumor/%s.hg19.%s_%s.masked_cnv.bed' % (line[1],line[0],line[0],min_cnv,min_size)],stderr=subprocess.STDOUT)
      if len(res.strip()) == 0:
        continue
    except subprocess.CalledProcessError as error:
      continue 
    output_dir = "output/TCGA/%s/%s/%s" % (line[1],line[0],'tumor')
    output.write("dsub --image us.gcr.io/aa-test-175718/aa --preemptible --project aa-test-175718 --zones \"us-west1-*\" --logging gs://aa-data-repo/%s/logging/ \\\n\t--output LOG_FILE=gs://aa-data-repo/%s/log \\\n\t --input-recursive DATA_REPO=gs://aa-data-repo/data_repo/ \\\n\t--output-recursive OUTPUT_DIR=gs://aa-data-repo/%s --disk-size 750 \\\n\t--input SH=%s BAM_FILE=%s BAI_FILE=%s.bai BED_FILE=gs://aa-data-repo/TCGA/%s/%s/tumor/%s.hg19.%s_%s.masked_cnv.bed \\\n\t--env PREFIX=%s --command 'sh ${SH}'\n" %  (output_dir,output_dir,output_dir,"gs://aa-data-repo/data_repo/test.sh",line[3],line[3],line[1],line[0],line[0],min_cnv,min_size,line[0]))
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

service = get_unauthorized_service(api='isb_cgc_tcga_api')
run_samples_python()
