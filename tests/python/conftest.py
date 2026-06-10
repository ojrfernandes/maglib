import sys
import os

# Make 'import maglib' work by adding the python/ source directory to sys.path.
# The compiled _maglib extension is placed in python/maglib/ by CMake.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

# Set non-interactive backend before any test imports matplotlib.pyplot.
# plot.py is imported at maglib package load time, so this must happen here.
import matplotlib
matplotlib.use('Agg')
