"""Top-level package for the SST-1M analysis and reconstruction pipeline."""

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _version

try:
    __version__ = _version(__name__)
except PackageNotFoundError:
    pass
