import wx
import wx.grid as grid
from wx.lib.plot import PolyLine, PlotCanvas, PlotGraphics
import matplotlib as plt
import numpy as np
plt.use('WXAgg')

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

from multipanel import *

class MyApp(wx.App):
    def __init__(self):
        super().__init__(clearSigInt=True)

        # init frame
        self.InitFrame()

    def InitFrame(self):
        frame = MainFrame(parent=None, title="My Frame", pos=(100, 100))
        frame.Show()


class MainFrame(wx.Frame):
    def __init__(self, parent, pos, title=wx.EmptyString, size = (800,600)):
        super(MainFrame, self).__init__(parent, size=(800,600), style=wx.DEFAULT_FRAME_STYLE)

        self.SetSizeHints( wx.DefaultSize, wx.DefaultSize )
        panelSizer = wx.GridSizer( 2, gap=(10, 10) )
        self.SetSizer( panelSizer )
        self.Layout()

        self.PanelOne = PanelOne(self)
        self.PanelTwo = PanelTwo(self)

        #self.contourPanel.Hide()
        #self.Centre( wx.BOTH )

    def __del__( self ):
        pass

    def loadNewPanel(self,invoker):
        if isinstance(invoker,PanelOne):
            #print "loading panel two"
            self.panel = PanelTwo(self)
        else:
            #print "loading panel one"
            self.panel = PanelOne(self)
            self.sizer_1.Fit(self)
            self.sizer_1.SetSizeHints(self)


class controlPanel(wx.Panel):
    def __init__(self, parent):
        super(controlPanel, self).__init__(parent)

        # domain table sizer
        self.domainSizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.domainSizer)

        # add text
        self.domainLabel = wx.StaticText(self, label = "Geometry Definition")
        #domainGrid.SetCellFont(0, 0, wx.Font(12, wx.ROMAN, wx.ITALIC, wx.NORMAL))
        self.domainLabel.SetFont(wx.Font(10, wx.SWISS, wx.NORMAL, wx.BOLD))
        self.domainSizer.Add(self.domainLabel)

        # geometry table
        self.domainGrid = grid.Grid(self)
        self.domainGrid.CreateGrid(6, 1)

        # set row and column values in domain table
        #domainGrid.SetCellFont(0, 0, wx.Font(12, wx.ROMAN, wx.ITALIC, wx.NORMAL))
        self.domainGrid.SetColLabelSize(0)
        self.domainGrid.SetCellValue(0, 0, "1.5")
        self.domainGrid.SetCellValue(1, 0, "1.3")
        self.domainGrid.SetCellValue(2, 0, "0.5")
        self.domainGrid.SetCellValue(3, 0, "20")
        self.domainGrid.SetCellValue(4, 0, "30")
        self.domainGrid.SetCellValue(5, 0, "26")

        self.domainGrid.SetColSize(0, 46)

        # set cell editors for M and N
        #domainGrid.SetCellEditor(4, 0, grid.GridCellNumberEditor(4, 0))
        #domainGrid.SetCellEditor(5, 0, grid.GridCellNumberEditor(5, 0))

        # set rownames
        self.domainGrid.SetRowLabelValue(0, "Length")
        self.domainGrid.SetRowLabelValue(1, "Height")
        self.domainGrid.SetRowLabelValue(2, "Wedge Start")
        self.domainGrid.SetRowLabelValue(3, "Wedge Angle")
        self.domainGrid.SetRowLabelValue(4, "M")
        self.domainGrid.SetRowLabelValue(5, "N")

        self.domainSizer.Add(self.domainGrid, wx.EXPAND, border=20)




    #     # add text
    #     sizer = wx.BoxSizer(wx.HORIZONTAL)
    #     self.label = wx.StaticText(self, label = "Hello")
    #     sizer.Add(self.label, 1, wx.EXPAND)


    #     # add button
    #     self.btn = wx.Button(self, label = "Click Here")
    #     self.btn.Bind(wx.EVT_BUTTON, self.onClick)
    #     sizer.Add(self.btn, 1)

    #     self.SetSizer(sizer)


    #     # add checkbox
    #     self.check1 = wx.CheckBox(self, label = 'Checc')
    #     self.Bind(wx.EVT_CHECKBOX, self.onCheck)
    #     sizer.Add(self.check1)


    # def onCheck(self, event):
    #     cb = event.GetEventObject()
    #     self.label.SetLabelText("Selected" + cb.GetLabel())
        
        

    # def onClick(self, event):
    #     self.label.SetLabelText("Text has been changed")


        # gridSizer = wx.GridSizer(4, 4, 5, 5)

        # for i in range(1, 17):
        #     btn = "Button" + str(i)

        #     gridSizer.Add(wx.Button(self,label=btn), 0, wx.EXPAND)
        #     self.SetSizer(gridSizer)

        # # box sizer
        # vbox = wx.BoxSizer(wx.VERTICAL)
        # hbox = wx.BoxSizer(wx.HORIZONTAL)

        # # add message to panel
        # text = wx.StaticText(self, id=wx.ID_ANY, label = 'yeeet', style = wx.ALIGN_CENTER_HORIZONTAL)
        # # ID_ANY means we don't care about the ID

        # vbox.Add(text, 0, wx.EXPAND)
        # self.SetSizer(vbox)

        # text2 = wx.StaticText(self, id=wx.ID_ANY, label = 'yeeet2', style = wx.ALIGN_CENTER_HORIZONTAL)
        # hbox.Add(text2, 0, wx.EXPAND)
        # vbox.Add(hbox)
        # self.SetSizer(vbox)

        # # add button here
        # button = wx.Button(parent=self, label="Clicc", pos = (20, 80))
        # button.Bind(event=wx.EVT_BUTTON, handler=self.onSubmit)


    # def onSubmit(self, event):
    #     # add button action
    #     webbrowser.open('https://google.com')


class contourPanel(wx.Panel): 
    def __init__(self, parent):
        super(contourPanel, self).__init__(parent)

            # create contour graph
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        #self.domainSizer.Add(self.canvas, 100, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()

        def draw(self):
            t = np.arange(0.0, 3.0, 0.01)
            s = np.sin(2 * np.pi * t)
            self.axes.plot(t, s)


if __name__== "__main__":
    app = MyApp()
    app.MainLoop()

