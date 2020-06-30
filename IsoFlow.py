import wx

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
## Start the application
###########################################################################
app.MainLoop()