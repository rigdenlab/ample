
import os

__all__ = ["AMPLE_DIR"]


AMPLE_DIR = os.sep.join(os.path.abspath(os.path.dirname(__file__)).split(os.sep)[ :-1 ])