import wx

import logging
import logging.config

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar

###########################################################################
## import wxFormBuilder file, initialize table values
###########################################################################	
import gui1

# initialize main frame  
app = wx.App(False, useBestVisual=True)
frame = gui1.MainFrame(None)
frame.SetTitle('IsoFlow 0.3.0')
#frame.SetIcon(wx.Icon("/Users/Elliot/Documents/Python/compressible_wedge/bmp/IsoFlow.ico"))
frame.SetPosition(wx.Point(240, 0))
frame.Show(True)

frame.contourPanel.figure = plt.figure( dpi=100, figsize=(5.6, 4.2), facecolor=(222/256,222/256,222/256) )
frame.contourPanel.cax = frame.contourPanel.figure.gca()

frame.iterPanel.figure = plt.figure( dpi=100, figsize=(5.1, 1.5), facecolor=(222/256,222/256,222/256) )
frame.iterPanel.iax = frame.iterPanel.figure.gca()


###########################################################################
## Start the application
###########################################################################
app.MainLoop()