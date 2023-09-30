import warnings
warnings.filterwarnings("default", category=DeprecationWarning)
def deprecated(msg):
    warnings.warn(msg, DeprecationWarning, stacklevel=2)
