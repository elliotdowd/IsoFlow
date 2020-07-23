import numpy as np
import matplotlib.pyplot as plt

# import domain values
M = 42
N = 18
length = 4
height = 1
theta1 = 60*(3.1415927/180)
y0 = 0.4
Aratio = 2.9
throat = 0.25

# Foesch nozzle parameters (see NAVAL ORDNANCE LABORATORY MEMORANDUM 10594 )

r0 = y0 / theta1
y1 = y0 * (np.sin(theta1)/theta1) * Aratio
x1 = (3/2) * (y1-y0) * (1/np.tan(theta1))

Nhalf = int(N/2)+1

# mesh left side of domain

x_i = np.linspace(0, (1+(1/M)), int(M/2)+3)
x_i = x_i * np.sin(x_i)**(2*throat/height)
y_i = np.linspace(-(1/2)*(1+(1/N)), (1/2)*(1+(1/N)), N+3)

xL = x_i / np.max(x_i)

y = y1 + (np.tan(theta1)/x1)* (xL**2) * (1-xL/(3*x1))
yL = ( y-np.min(y) )
yL = ( yL / np.max(yL) )

xx1, yy1 = np.meshgrid(xL*length, y_i*height)
xx1 = np.transpose(xx1)
yy1 = np.transpose(yy1)

yy1 = np.fliplr(yy1)

for j in range(0, Nhalf+1):
    yy1[:,Nhalf-j] = yy1[:,Nhalf-j] + (1/2)*((height/throat)-1)*yL*((j)/Nhalf)
    yy1[:,Nhalf+j] = yy1[:,Nhalf+j] - (1/2)*((height/throat)-1)*yL*((j)/Nhalf)

yy1 = (height / np.max(yy1[:,1:-1])) * yy1

xx1 = np.vstack( [np.flipud(-xx1[-1,:]-(1/M)*length), np.flipud(-xx1), xx1[1:,:]] )
yy1 = np.vstack( [(yy1[-1,:]), np.flipud(yy1), yy1[1:,:]] )

yy1 = np.fliplr(yy1)

# plt.plot(xL, yL)
plt.plot(xx1, yy1, 'k')
plt.plot(np.transpose(xx1), np.transpose(yy1), 'k')

plt.show()