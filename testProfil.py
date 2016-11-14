import numpy as np
import triangleLib as tl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import time
import sys
sys.path.append('./VortexLatticeMethod/')
sys.path.append('./asciiArt/')
sys.path.append('./profiles/')
import splineProfileMultiParam as prf
import vlm
import asciiArt
import sys
import mesh




asciiArt.openFile('./asciiArt/minusEtCortex')

pf0 = prf.Profile(typ = 'fon',par = [0.82,0.21,0.13,0.04,0.029],npt = 5)
mesh.thickWing(pf0,'wing2D')
asciiArt.openFile('./asciiArt/snoopy')
