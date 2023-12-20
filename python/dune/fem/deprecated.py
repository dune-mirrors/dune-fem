import warnings
"""
warnings.simplefilter("always", category=DeprecationWarning)

# ignore a specific warning caused by the plotting functions:
# matplotlib_inline/config.py:68: DeprecationWarning: InlineBackend._figure_format_changed is deprecated in traitlets 4.1: use @observe and @unobserve instead.
#   def _figure_format_changed(self, name, old, new):
warnings.filterwarnings(
    action='ignore',
    category=DeprecationWarning,
    module=r'.*matplotlib_inline/config'
)
"""
warnings.filterwarnings(
    action='always',
    category=DeprecationWarning,
    module=r'dune'
)

def deprecated(msg):
    msg += """
To disable this message add 'warnings.filterwarnings("ignore", category=DeprecationWarning)'
to your script after importing all dune.fem modules.
"""
    warnings.warn(msg, DeprecationWarning, stacklevel=3)
