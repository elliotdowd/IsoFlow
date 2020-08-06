import numpy as np
from stl import mesh

import matplotlib.pyplot as plt

square = mesh.Mesh.from_file('python\mesh\grid\square.stl')

square.normals