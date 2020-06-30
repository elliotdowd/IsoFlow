import wx

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar

import logging
import logging.config

###########################################################################
## import wxFormBuilder file, initialize table values
###########################################################################	
import gui1

# initialize main frame  
app = wx.App(False, useBestVisual=True)
frame = gui1.MainFrame(None)
frame.SetTitle('Flux Vector Differencing CFD Solver')
frame.SetPosition(wx.Point(240, 0))
frame.Show(True)

###########################################################################
## Plotting Panels
###########################################################################	

frame.contourPanel.figure = plt.figure( dpi=100, figsize=(5.5, 3.8), facecolor=(222/256,222/256,222/256) )
frame.contourPanel.cax = frame.contourPanel.figure.gca()
frame.contourPanel.cax.set_position([0.1, 0.18, 0.84, 0.82])

frame.iterPanel.figure = plt.figure( dpi=100, figsize=(5, 1), facecolor=(222/256,222/256,222/256) )
frame.iterPanel.iax = frame.iterPanel.figure.gca()

###########################################################################
## Start the application
###########################################################################
app.MainLoop()