import pkg_resources

try:
    __VERSION__ = pkg_resources.get_distribution('pbcoretools').version
except Exception:
    __VERSION__ = 'unknown'
