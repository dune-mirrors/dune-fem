import warnings
# warnings.filterwarnings("default", category=DeprecationWarning)
# warnings.simplefilter("default", category=DeprecationWarning)
warnings.simplefilter("always", category=DeprecationWarning)
def deprecated(msg):
    msg += """
To disable this message add 'warnings.filterwarnings("ignore", category=DeprecationWarning)'
to your script after importing all dune.fem modules.
"""
    warnings.warn(msg, DeprecationWarning, stacklevel=3)
