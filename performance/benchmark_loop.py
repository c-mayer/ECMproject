#!/usr/bin/env python

"""Python script to loop benchmark_test.py

@author: Marcus Holzer, Christian Mayer
"""

import subprocess

model_numbers = range(16,18)
repeats = [3,3] #* len(model_numbers) # give list of repeats with same length as model_numbers
n_cores = 20
tools = 'own'
for i, r in zip(model_numbers, repeats):
    cmd = ['./benchmark_test.py', '-f', '../metabolic_models/mmsyn_sm' + str(i).zfill(2) + '.xml', '-n', str(n_cores), '-r', str(r), '-m', 'mmsyn_sm' + str(i).zfill(2), '-t', tools, '-op', './times/ECMproject/single/']
    subprocess.run(cmd)
