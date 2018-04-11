# Script for multiple runs of wanderlust 
import warnings

# data processing
with warnings.catch_warnings():
    warnings.simplefilter('ignore', FutureWarning)
    import scanpy.api as sc
import numpy as np
import scanpy.api as sc
import pandas as pd
import os
import scipy.sparse as sp

# plotting
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


sns.set_style('white')
matplotlib.rcParams['figure.figsize'] = [4, 4]
matplotlib.rcParams['figure.dpi'] = 100
warnings.filterwarnings(action="ignore", module="matplotlib", message="findfont")

# change logging settings for scanpy
sc.settings.verbosity = 4  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80)  # low dpi (dots per inch) yields small inline figures
sc.settings.n_jobs = 30
sc.logging.print_version_and_date()
sc.logging.print_versions_dependencies_numerics()