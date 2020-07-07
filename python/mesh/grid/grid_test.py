import numpy as np
import matplotlib.pyplot as plt

# import domain values
M = 24
N = 84
length = 3
height = 1
theta1 = 0
cyl_start = 1
cyl_end = 2


r = np.linspace( cyl_start*(1-(1/M)), length*(1+(1/M)), M+3)
th = np.linspace( 0, 2*np.pi, N+3)


rr, tt = np.meshgrid(r, th)
xx = rr * np.cos(tt)
yy = rr * np.sin(tt)

xx = np.transpose(xx)
yy = np.transpose(yy)

fig = plt.figure( )
ax = plt.subplot(111)

mesh = ax.plot( xx, yy, 'r-'); ax.plot( np.transpose(xx), np.transpose(yy), 'r-' )

ax.set_aspect('equal')
plt.show()