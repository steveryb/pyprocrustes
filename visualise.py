import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import Procrustes
def load_curve_data(filename):
    with open(filename) as f:
        line = f.readline()
        data = line.split(",")
        version = data[0]
        num_points = data[1]
        data = data[2:]
        points = [tuple((float(data[i+j]) for j in range(3))) for i in range(0,len(data),3)]
        return points
    

def draw_curves(curves_to_draw):
    mpl.rcParams['legend.fontsize'] = 10
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for curve in curves_to_draw:
        x,y,z = [[point[i] for point in curve] for i in range(3)]
        ax.plot(x,y,zs=z, label="curve")
        ax.legend()
    plt.show()

def visualise(filenames):
    curves_to_draw = [load_curve_data(filename) for filename in filenames]
    draw_curves(curves_to_draw)

curves = [load_curve_data(f) for f in ["curves/13 1","curves/ref_1.crv"]]
curves = Procrustes.superposition(*[np.array(c) for c in curves])
print("distance",Procrustes.min_distance(*curves))
draw_curves(curves)
