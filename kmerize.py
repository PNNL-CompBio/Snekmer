#!/usr/bin/env python

import os, sys, gzip, random, copy
from types import *

from Util.Options import *
from Util.SIEVEInit import *

    def reduce_alphabet(character, mapping=None, **kw):
        for key in mapping.keys():
            if character in key:
                return mapping[key]
