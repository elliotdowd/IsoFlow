# -*- coding: utf-8 -*- 

###########################################################################
## Python code generated with wxFormBuilder (version Jun 17 2015)
## http://www.wxformbuilder.org/
##
## PLEASE DO "NOT" EDIT THIS FILE!
###########################################################################

import wx
import wx.xrc
import wx.grid

class MyApp(wx.App):
    def __init__(self):
        super().__init__(clearSigInt=True)

        # init frame
        self.InitFrame()

    def InitFrame(self):
        frame = MainFrame(parent=None, pos=(100, 100))
        frame.Show()

###########################################################################
## Class MainFrame
###########################################################################

class MainFrame ( wx.Frame ):
	
	def __init__( self, parent ):
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = wx.EmptyString, pos = wx.DefaultPosition, size = wx.Size( 800,900 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
		
		self.SetSizeHintsSz( wx.DefaultSize, wx.DefaultSize )
		
		MainSizer = wx.GridBagSizer( 0, 0 )
		MainSizer.SetFlexibleDirection( wx.BOTH )
		MainSizer.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
		
		self.domainGrid = wx.grid.Grid( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, 0 )
		
		# Grid
		self.domainGrid.CreateGrid( 6, 1 )
		self.domainGrid.EnableEditing( True )
		self.domainGrid.EnableGridLines( True )
		self.domainGrid.EnableDragGridSize( False )
		self.domainGrid.SetMargins( 0, 0 )
		
		# Columns
		self.domainGrid.EnableDragColMove( False )
		self.domainGrid.EnableDragColSize( True )
		self.domainGrid.SetColLabelSize( 0 )
		self.domainGrid.SetColLabelAlignment( wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
		
		# Rows
		self.domainGrid.EnableDragRowSize( True )
		self.domainGrid.SetRowLabelSize( 102 )
		self.domainGrid.SetRowLabelValue( 0, u"Length" )
		self.domainGrid.SetRowLabelValue( 1, u"Height" )
		self.domainGrid.SetRowLabelValue( 2, u"Wedge Start" )
		self.domainGrid.SetRowLabelValue( 3, u"Wedge Angle" )
		self.domainGrid.SetRowLabelValue( 4, u"Horiz. Cells" )
		self.domainGrid.SetRowLabelValue( 5, u"Vertical Cells" )
		self.domainGrid.SetRowLabelAlignment( wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
		
		# Label Appearance
		
		# Cell Defaults
		self.domainGrid.SetDefaultCellAlignment( wx.ALIGN_LEFT, wx.ALIGN_TOP )
		MainSizer.Add( self.domainGrid, wx.GBPosition( 1, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		self.m_staticText1 = wx.StaticText( self, wx.ID_ANY, u"Geometry Parameters", wx.DefaultPosition, wx.DefaultSize, 0 )
		self.m_staticText1.Wrap( -1 )
		self.m_staticText1.SetFont( wx.Font( 10, 74, 90, 92, False, "Arial" ) )
		
		MainSizer.Add( self.m_staticText1, wx.GBPosition( 0, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		self.contourPanel = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		MainSizer.Add( self.contourPanel, wx.GBPosition( 1, 2 ), wx.GBSpan( 1, 3 ), wx.EXPAND |wx.ALL, 5 )
		
		self.m_button1 = wx.Button( self, wx.ID_ANY, u"Generate Grid", wx.DefaultPosition, wx.DefaultSize, 0 )
		MainSizer.Add( self.m_button1, wx.GBPosition( 3, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		meshTypeChoices = [ u"Wedge", u"Airfoil" ]
		self.meshType = wx.ListBox( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, meshTypeChoices, 0 )
		MainSizer.Add( self.meshType, wx.GBPosition( 2, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		
		self.SetSizer( MainSizer )
		self.Layout()
		
		self.Centre( wx.BOTH )
	
	def __del__( self ):
		pass
	

if __name__== "__main__":
    app = MyApp()
    app.MainLoop()
