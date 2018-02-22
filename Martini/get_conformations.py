#!/sansom/s137/bioc1535/anaconda2/bin/python

from __future__ import print_function, division

import os
import re
import shlex
import sys
import gromacs
import numpy as np

def parse_xvg(fname, sel_columns='all'):
    """Parses XVG file legends and data"""
    
    _ignored = set(('legend', 'view'))
    _re_series = re.compile('s[0-9]+$')
    _re_xyaxis = re.compile('[xy]axis$')

    metadata = {}
    num_data = []
    
    metadata['labels'] = {}
    metadata['labels']['series'] = []

    ff_path = os.path.abspath(fname)
    if not os.path.isfile(ff_path):
        raise IOError('File not readable: {0}'.format(ff_path))
    
    with open(ff_path, 'r') as fhandle:
        for line in fhandle:
            line = line.strip()
            if line.startswith('@'):
                tokens = shlex.split(line[1:])
                if tokens[0] in _ignored:
                    continue
                elif tokens[0] == 'TYPE':
                    if tokens[1] != 'xy':
                        raise ValueError('Chart type unsupported: \'{0}\'. Must be \'xy\''.format(tokens[1]))
                elif _re_series.match(tokens[0]):
                    metadata['labels']['series'].append(tokens[-1])
                elif _re_xyaxis.match(tokens[0]):
                    metadata['labels'][tokens[0]] = tokens[-1]
                elif len(tokens) == 2:
                    metadata[tokens[0]] = tokens[1]
                else:
                    print('Unsupported entry: {0} - ignoring'.format(tokens[0]), file=sys.stderr)
            elif line[0].isdigit():
                numbers = line.split()
                num_data.append((numbers[0], numbers[1]))
    
    if sel_columns != 'all':
        sel_columns = map(int, sel_columns)
        x_axis = num_data[0]
        num_data = [x_axis] + [num_data[col] for col in sel_columns]
        metadata['labels']['series'] = [metadata['labels']['series'][col - 1] for col in sel_columns]
    
    return metadata, num_data

if __name__ == '__main__':
	
    import argparse
    from argparse import RawDescriptionHelpFormatter
    
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)

    ap.add_argument('xvg_f', type=str, help='XVG input file', metavar='XVG input file')

    cmd = ap.parse_args()
    metadata, data = parse_xvg(cmd.xvg_f)

    init_val = 1.5 #CHANGE AS APPROPRIATE - Should be a multiple of 0.1
    end_val = 1.7 #CHANGE AS APPROPRIATE - Should be a multiple of 0.1

    num_windows = int(10*(end_val-init_val))
    difference = 0.05	
	
    window_id = 1
    overall_dict = {}	
    current_best_window = 0
    current_best_value = 0
	
    for j in range(0,num_windows+1):
        for i in range(0,len(data)):
            if abs(float(data[i][1]) - init_val) < difference:
                current_best_window = float(data[i][0])
                current_best_value = float(data[i][1])
                difference = abs(float(data[i][1]) - init_val)
				
        overall_dict[window_id] = current_best_window
        cmd = gromacs.trjconv(f='pull.xtc', o='conf'+str(window_id)+'.pdb', b=str(current_best_window), e=str(current_best_window), s='pull.tpr',  input=('system'))
        window_id += 1
        init_val += 0.1
        difference = 0.05
        current_best_window = 0
        current_best_value = 0	
