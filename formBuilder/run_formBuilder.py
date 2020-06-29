import wx

import matplotlib as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar
  
# import the newly created GUI file 
import gui1
		
        
app = wx.App(False)
frame = gui1.MainFrame(None)
frame.SetPosition(wx.Point(240, 0))
frame.Show(True)

###########################################################################
## Plotting Panels
###########################################################################	
frame.contourPanel.figure = plt.figure.Figure( dpi=None, figsize=(2, 2) )
frame.contourPanel.cax = frame.contourPanel.figure.add_subplot()
frame.contourPanel.cax.plot( [1, 2, 3], [2, 1, 4] )
frame.contourPanel.canvas = FigureCanvas(frame.contourPanel, -1, frame.contourPanel.figure)


#start the application
app.MainLoop() 