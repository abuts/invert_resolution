import os
from mantid.simpleapi import *
from mantid import api
import numpy as np
import tube
ws = Load("MER31013.n003")
ws = Integration(ws,1500,9000)
pos = [0.84,0.44,0.02,-0.45,-0.65]
tt = [1,1,1,1,1]
calt = tube.calibrate(ws,"MERLIN",pos,tt)
