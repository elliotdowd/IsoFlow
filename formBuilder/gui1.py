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

###########################################################################
## Class MainFrame
###########################################################################

class MainFrame ( wx.Frame ):
	
	def __init__( self, parent ):
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = wx.EmptyString, pos = wx.Point( 100,100 ), size = wx.Size( 800,900 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
		
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
		self.domainGrid.SetColSize( 0, 72 )
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
		MainSizer.Add( self.domainGrid, wx.GBPosition( 2, 0 ), wx.GBSpan( 1, 1 ), wx.ALL|wx.EXPAND, 5 )
		
		self.contourPanel = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		MainSizer.Add( self.contourPanel, wx.GBPosition( 1, 2 ), wx.GBSpan( 7, 54 ), wx.ALL|wx.EXPAND, 5 )
		
		self.iterPanel = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		MainSizer.Add( self.iterPanel, wx.GBPosition( 11, 2 ), wx.GBSpan( 3, 54 ), wx.ALL|wx.EXPAND, 5 )
		
		self.parameterGrid = wx.grid.Grid( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, 0 )
		
		# Grid
		self.parameterGrid.CreateGrid( 3, 1 )
		self.parameterGrid.EnableEditing( True )
		self.parameterGrid.EnableGridLines( True )
		self.parameterGrid.EnableDragGridSize( False )
		self.parameterGrid.SetMargins( 0, 0 )
		
		# Columns
		self.parameterGrid.SetColSize( 0, 54 )
		self.parameterGrid.EnableDragColMove( False )
		self.parameterGrid.EnableDragColSize( True )
		self.parameterGrid.SetColLabelSize( 0 )
		self.parameterGrid.SetColLabelAlignment( wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
		
		# Rows
		self.parameterGrid.EnableDragRowSize( True )
		self.parameterGrid.SetRowLabelSize( 120 )
		self.parameterGrid.SetRowLabelValue( 0, u"Inlet Mach #" )
		self.parameterGrid.SetRowLabelValue( 1, u"Inlet Pres. (Pa)" )
		self.parameterGrid.SetRowLabelValue( 2, u"Inlet Temp. (K)" )
		self.parameterGrid.SetRowLabelAlignment( wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
		
		# Label Appearance
		
		# Cell Defaults
		self.parameterGrid.SetDefaultCellAlignment( wx.ALIGN_LEFT, wx.ALIGN_TOP )
		MainSizer.Add( self.parameterGrid, wx.GBPosition( 7, 0 ), wx.GBSpan( 1, 1 ), wx.ALL|wx.EXPAND, 5 )
		
		self.simGrid = wx.grid.Grid( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, 0 )
		
		# Grid
		self.simGrid.CreateGrid( 3, 1 )
		self.simGrid.EnableEditing( True )
		self.simGrid.EnableGridLines( True )
		self.simGrid.EnableDragGridSize( False )
		self.simGrid.SetMargins( 0, 0 )
		
		# Columns
		self.simGrid.SetColSize( 0, 54 )
		self.simGrid.EnableDragColMove( False )
		self.simGrid.EnableDragColSize( True )
		self.simGrid.SetColLabelSize( 0 )
		self.simGrid.SetColLabelAlignment( wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
		
		# Rows
		self.simGrid.EnableDragRowSize( True )
		self.simGrid.SetRowLabelSize( 120 )
		self.simGrid.SetRowLabelValue( 0, u"Max. CFL" )
		self.simGrid.SetRowLabelValue( 1, u"Inlet Pres. (Pa)" )
		self.simGrid.SetRowLabelValue( 2, u"Residual Tol." )
		self.simGrid.SetRowLabelValue( 3, u"Iterations" )
		self.simGrid.SetRowLabelValue( 4, wx.EmptyString )
		self.simGrid.SetRowLabelAlignment( wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
		
		# Label Appearance
		
		# Cell Defaults
		self.simGrid.SetDefaultCellAlignment( wx.ALIGN_LEFT, wx.ALIGN_TOP )
		MainSizer.Add( self.simGrid, wx.GBPosition( 11, 0 ), wx.GBSpan( 1, 1 ), wx.ALL|wx.EXPAND, 5 )
		
		self.m_staticText11 = wx.StaticText( self, wx.ID_ANY, u"Inlet Parameters", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_CENTRE )
		self.m_staticText11.Wrap( -1 )
		self.m_staticText11.SetFont( wx.Font( 10, 74, 90, 92, False, "Arial" ) )
		
		MainSizer.Add( self.m_staticText11, wx.GBPosition( 6, 0 ), wx.GBSpan( 1, 1 ), wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.m_staticText111 = wx.StaticText( self, wx.ID_ANY, u"Simulation Options", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_CENTRE )
		self.m_staticText111.Wrap( -1 )
		self.m_staticText111.SetFont( wx.Font( 10, 74, 90, 92, False, "Arial" ) )
		
		MainSizer.Add( self.m_staticText111, wx.GBPosition( 10, 0 ), wx.GBSpan( 1, 1 ), wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		gridChoiceChoices = [ u"Wedge", u"Airfoil" ]
		self.gridChoice = wx.Choice( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, gridChoiceChoices, 0 )
		self.gridChoice.SetSelection( 0 )
		MainSizer.Add( self.gridChoice, wx.GBPosition( 3, 0 ), wx.GBSpan( 1, 1 ), wx.ALL|wx.EXPAND, 5 )
		
		self.m_button3 = wx.Button( self, wx.ID_ANY, u"Run Simulation", wx.DefaultPosition, wx.DefaultSize, 0 )
		MainSizer.Add( self.m_button3, wx.GBPosition( 13, 0 ), wx.GBSpan( 1, 1 ), wx.ALL|wx.EXPAND, 5 )
		
		schemeChoiceChoices = [ u"AUSM", u"AUSM+up", u"AUSMDV" ]
		self.schemeChoice = wx.Choice( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, schemeChoiceChoices, 0 )
		self.schemeChoice.SetSelection( 0 )
		MainSizer.Add( self.schemeChoice, wx.GBPosition( 12, 0 ), wx.GBSpan( 1, 1 ), wx.ALL|wx.EXPAND, 5 )
		
		self.gridButton = wx.Button( self, wx.ID_ANY, u"Generate Grid", wx.DefaultPosition, wx.DefaultSize, 0 )
		MainSizer.Add( self.gridButton, wx.GBPosition( 4, 0 ), wx.GBSpan( 1, 1 ), wx.ALL|wx.EXPAND, 5 )
		
		self.initButton = wx.Button( self, wx.ID_ANY, u"Initialize", wx.DefaultPosition, wx.DefaultSize, 0 )
		MainSizer.Add( self.initButton, wx.GBPosition( 8, 0 ), wx.GBSpan( 1, 1 ), wx.ALL|wx.EXPAND, 5 )
		
		self.m_staticText1 = wx.StaticText( self, wx.ID_ANY, u"Geometry Parameters", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_CENTRE )
		self.m_staticText1.Wrap( -1 )
		self.m_staticText1.SetFont( wx.Font( 10, 74, 90, 92, False, "Arial" ) )
		
		MainSizer.Add( self.m_staticText1, wx.GBPosition( 1, 0 ), wx.GBSpan( 1, 1 ), wx.ALL, 5 )
		
		
		self.SetSizer( MainSizer )
		self.Layout()
		
		self.Centre( wx.BOTH )
	
	def __del__( self ):
		pass
	

