#!/usr/bin/env python3

import pandas as pd
import os
import psutil

d = pd.read_csv( './SDSS_SQL_data_out.txt', skiprows=1 )
#print( d['ra'][:10] )

p = psutil.Process( os.getpid() )
x = p.memory_info().rss / 1024 / 1024
print( x )
