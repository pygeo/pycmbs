#!/bin/bash
#preprocess all data

python preprocess_cmip5_offline.py amip rsds
python preprocess_cmip5_offline.py amip rsus
python preprocess_cmip5_offline.py historical rsds
python preprocess_cmip5_offline.py historical rsus

