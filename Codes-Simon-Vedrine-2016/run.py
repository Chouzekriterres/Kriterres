#!/usr/bin/python3

import sys
import os
import argparse

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import utils as u
#import processing_tools as tools
from processing_tools import *

dpath = os.path.abspath('data') + '/'

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('filename', help='name of data file including path')
	args = parser.parse_args()
	
	filename = args.filename
	
	if '.ASC' in filename:
		data = u.read_ascii('%s%s' % (dpath, filename))
		dx_dt = (0.00667, 0.0185)
	
	elif '.out' in filename:
		data, dx_dt = u.read_hdf5('%s%s' % (dpath, filename))
		user_gain (data, dx_dt, 'linear', (2,0), (0,100), plot=True)
                
	else:
		data, dx_dt = u.read_ramac('%s%s' % (dpath, filename))
		
		tools.save_new_version(data, 'allo.bin')
	
	
if __name__ == '__main__':
	main()
