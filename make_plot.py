from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import ticker
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

import numpy

ticker.rcParams['xtick.direction'] = 'in'
ticker.rcParams['ytick.direction'] = 'in'

lines = []

f = open("out.txt")

for line in f.readlines():
	lines.append(line)

X = numpy.linspace(0, 1, num=len(lines))
Y = numpy.linspace(0, 1, num=len(lines))

X, Y = numpy.meshgrid(X, Y)

Z = numpy.array([[0 for i in range(len(lines))] for j in range(len(lines))])

for i in range(len(lines)):
	line = lines[i]

	f_line = [float(x) for x in line.split(" ")[:-2]]

	for j in range(len(f_line)):
		Z[i][j] = f_line[j]

f.close()

fig1 = Figure()
FigureCanvas(fig1)

axs1 = fig1.add_subplot(111, projection='3d')

axs1.plot_surface(X, Y, Z, linewidth=0, antialiased=False, cmap=cm.coolwarm)
axs1.set_xlim([0, 1])
axs1.set_ylim([0, 1])

fig1.savefig("plot.png", fmt="png")