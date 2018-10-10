#!/usr/bin/env python
# Make sure you already loaded the 'Python' module (python 3)

import sj2psi as sj
import tempfile

def get_usage(filename):
  # Convert this to a Dict so it can be imported into R code via rPython
  df = sj.get_psis(sj.read_sj_out_tab(filename), min_unique=0, min_multimap=0)
  fh, fn = tempfile.mkstemp()
  df.to_csv(fn, index=False)
  return fn
