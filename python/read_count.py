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


PANCANCER_DIR=os.environ['PANCANCER']
TEMP_DIR= os.environ['TEMP_DIR'] if 'TEMP_DIR' in os.environ else None

