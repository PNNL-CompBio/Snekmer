import importlib.metadata

from . import alphabet
from . import utils
from . import io
from . import plot
from . import score
from . import model
from . import cluster
from . import vectorize
from . import report

__version__ = importlib.metadata.version("snekmer")
