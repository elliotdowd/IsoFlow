import wx

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar

###########################################################################
## import wxFormBuilder file, initialize table values
###########################################################################	
import gui1
		
# initialize main frame  
app = wx.App(False)
frame = gui1.MainFrame(None)
frame.SetTitle('Flux Vector Differencing CFD Solver')
frame.SetPosition(wx.Point(240, 0))
frame.Show(True)

# set domain row values
frame.domainGrid.SetCellValue( 0, 0, "1.5")
frame.domainGrid.SetCellValue( 1, 0, "1.3")
frame.domainGrid.SetCellValue( 2, 0, "0.5")
frame.domainGrid.SetCellValue( 3, 0, "1.1")
frame.domainGrid.SetCellValue( 4, 0, "20")
frame.domainGrid.SetCellValue( 5, 0, "36")
frame.domainGrid.SetCellValue( 6, 0, "32")

# set initialization row values
frame.parameterGrid.SetCellValue( 0, 0, "3.0")
frame.parameterGrid.SetCellValue( 1, 0, "101325")
frame.parameterGrid.SetCellValue( 2, 0, "300")

# set initialization row values
frame.simGrid.SetCellValue( 0, 0, "0.4")
frame.simGrid.SetCellValue( 1, 0, "1000")
frame.simGrid.SetCellValue( 2, 0, "-6")

###########################################################################
## Plotting Panels
###########################################################################	

# contour panel plotting
frame.contourPanel.figure = plt.figure( dpi=100, figsize=(5.4, 3.2) )
frame.contourPanel.cax = frame.contourPanel.figure.gca(projection='3d')
frame.contourPanel.cax.plot( [1, 2, 3], [2, 1, 4] )
frame.contourPanel.cax.set_xlabel('x-coordinate (m)')
frame.contourPanel.cax.set_ylabel('y-coordinate (m)')
frame.contourPanel.canvas = FigureCanvas(frame.contourPanel, -1, frame.contourPanel.figure)


frame.iterPanel.figure = mpl.figure.Figure( dpi=100, figsize=(5.4, 0.9) )
frame.iterPanel.iax = frame.iterPanel.figure.add_subplot(111)
frame.iterPanel.iax.plot( [1, 2, 3], [0, -1, -3] )
frame.contourPanel.canvas = FigureCanvas(frame.iterPanel, -1, frame.iterPanel.figure)


###########################################################################
## Button binding events
###########################################################################

def generateClick(self, event):
   btn = event.GetEventObject().GetLabel()
   print("Label of pressed button = "), btn

###########################################################################
## Start the application
###########################################################################
app.MainLoop()