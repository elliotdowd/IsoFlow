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
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = wx.EmptyString, pos = wx.Point( 100,100 ), size = wx.Size( 800,850 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
		
		self.SetSizeHints( wx.DefaultSize, wx.DefaultSize )
		self.SetForegroundColour( wx.SystemSettings.GetColour( wx.SYS_COLOUR_INACTIVEBORDER ) )
		self.SetBackgroundColour( wx.Colour( 222, 222, 222 ) )
		
		MainSizer = wx.GridBagSizer( 0, 0 )
		MainSizer.SetFlexibleDirection( wx.BOTH )
		MainSizer.SetNonFlexibleGrowMode( wx.FLEX_GROWMODE_SPECIFIED )
		
		self.domainGrid = wx.grid.Grid( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, 0 )
		
		# Grid
		self.domainGrid.CreateGrid( 7, 1 )
		self.domainGrid.EnableEditing( True )
		self.domainGrid.EnableGridLines( True )
		self.domainGrid.EnableDragGridSize( False )
		self.domainGrid.SetMargins( 0, 0 )
		
		# Columns
		self.domainGrid.SetColSize( 0, 32 )
		self.domainGrid.EnableDragColMove( False )
		self.domainGrid.EnableDragColSize( True )
		self.domainGrid.SetColLabelSize( 0 )
		self.domainGrid.SetColLabelAlignment( wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
		
		# Rows
		self.domainGrid.EnableDragRowSize( True )
		self.domainGrid.SetRowLabelSize( 112 )
		self.domainGrid.SetRowLabelValue( 0, u"Length" )
		self.domainGrid.SetRowLabelValue( 1, u"Height" )
		self.domainGrid.SetRowLabelValue( 2, u"Wedge Start" )
		self.domainGrid.SetRowLabelValue( 3, u"Wedge End" )
		self.domainGrid.SetRowLabelValue( 4, u"Wedge Angle" )
		self.domainGrid.SetRowLabelValue( 5, u"Horizontal Cells" )
		self.domainGrid.SetRowLabelValue( 6, u"Vertical Cells" )
		self.domainGrid.SetRowLabelAlignment( wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
		
		# Label Appearance
		
		# Cell Defaults
		self.domainGrid.SetDefaultCellAlignment( wx.ALIGN_LEFT, wx.ALIGN_TOP )
		MainSizer.Add( self.domainGrid, wx.GBPosition( 2, 0 ), wx.GBSpan( 1, 1 ), wx.ALL|wx.EXPAND, 5 )
		
		self.contourPanel = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.contourPanel.SetBackgroundColour( wx.Colour( 222, 222, 222 ) )
		
		MainSizer.Add( self.contourPanel, wx.GBPosition( 1, 2 ), wx.GBSpan( 8, 68 ), wx.ALL|wx.EXPAND, 5 )
		
		self.m_toolBar1 = wx.ToolBar( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TB_HORIZONTAL ) 
		self.m_toolBar1.Realize() 
		
		MainSizer.Add( self.m_toolBar1, wx.GBPosition( 0, 0 ), wx.GBSpan( 1, 62 ), wx.EXPAND, 5 )
		
		self.iterPanel = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.iterPanel.SetBackgroundColour( wx.Colour( 255, 255, 255 ) )
		
		MainSizer.Add( self.iterPanel, wx.GBPosition( 10, 2 ), wx.GBSpan( 4, 54 ), wx.ALL|wx.EXPAND, 5 )
		
		self.parameterGrid = wx.grid.Grid( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, 0 )
		
		# Grid
		self.parameterGrid.CreateGrid( 3, 1 )
		self.parameterGrid.EnableEditing( True )
		self.parameterGrid.EnableGridLines( True )
		self.parameterGrid.EnableDragGridSize( False )
		self.parameterGrid.SetMargins( 0, 0 )
		
		# Columns
		self.parameterGrid.SetColSize( 0, 44 )
		self.parameterGrid.EnableDragColMove( False )
		self.parameterGrid.EnableDragColSize( True )
		self.parameterGrid.SetColLabelSize( 0 )
		self.parameterGrid.SetColLabelAlignment( wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
		
		# Rows
		self.parameterGrid.EnableDragRowSize( True )
		self.parameterGrid.SetRowLabelSize( 100 )
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
		self.simGrid.SetColSize( 0, 44 )
		self.simGrid.EnableDragColMove( False )
		self.simGrid.EnableDragColSize( True )
		self.simGrid.SetColLabelSize( 0 )
		self.simGrid.SetColLabelAlignment( wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
		
		# Rows
		self.simGrid.EnableDragRowSize( True )
		self.simGrid.SetRowLabelSize( 100 )
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
		self.m_staticText11.SetForegroundColour(wx.Colour(0, 0, 0))
		
		MainSizer.Add( self.m_staticText11, wx.GBPosition( 6, 0 ), wx.GBSpan( 1, 1 ), wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.m_staticText111 = wx.StaticText( self, wx.ID_ANY, u"Simulation Options", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_CENTRE )
		self.m_staticText111.Wrap( -1 )
		self.m_staticText111.SetFont( wx.Font( 10, 74, 90, 92, False, "Arial" ) )
		self.m_staticText111.SetForegroundColour(wx.Colour(0, 0, 0))

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
		
		self.m_staticText1 = wx.StaticText( self, wx.ID_ANY, u"Domain Parameters", wx.DefaultPosition, wx.DefaultSize, wx.ALIGN_CENTRE )
		self.m_staticText1.Wrap( -1 )
		self.m_staticText1.SetFont( wx.Font( 10, 74, 90, 92, False, "Arial" ) )
		self.m_staticText1.SetForegroundColour(wx.Colour(0, 0, 0))
		
		MainSizer.Add( self.m_staticText1, wx.GBPosition( 1, 0 ), wx.GBSpan( 1, 1 ), wx.ALL|wx.ALIGN_CENTER_HORIZONTAL, 5 )
		
		self.m_staticline1 = wx.StaticLine( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LI_HORIZONTAL )
		MainSizer.Add( self.m_staticline1, wx.GBPosition( 5, 0 ), wx.GBSpan( 1, 1 ), wx.EXPAND |wx.ALL, 5 )
		
		self.m_staticline11 = wx.StaticLine( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.LI_HORIZONTAL )
		MainSizer.Add( self.m_staticline11, wx.GBPosition( 9, 0 ), wx.GBSpan( 1, 1 ), wx.EXPAND |wx.ALL, 5 )
		
		
		self.SetSizer( MainSizer )
		self.Layout()
		
		self.Centre( wx.BOTH )
		
		# Connect Events
		self.m_button3.Bind( wx.EVT_BUTTON, self.call_scheme )
		self.gridButton.Bind( wx.EVT_BUTTON, self.call_grid )
		self.initButton.Bind( wx.EVT_BUTTON, self.call_init )
	
	def __del__( self ):
		pass
	
	
	# Virtual event handlers, overide them in your derived class
	def call_scheme( self, event ):

		import numpy as np
		from pytictoc import TicToc
		from python.finite_volume.AUSM.schemes import AUSM, AUSMplusup, AUSMDV
		import matplotlib.pyplot as plt
		from matplotlib import cm
		import matplotlib as mpl
		from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

		t = TicToc()

		# run AUSM family scheme
		t.tic()
		scheme = self.schemeChoice.Strings[self.schemeChoice.Selection]
		
		if scheme == 'AUSM':
			self.state = AUSM( self.domain, self.mesh, self.parameters, self.state, self.gas )
		elif scheme == 'AUSM+up':
			self.state = AUSMplusup( self.domain, self.mesh, self.parameters, self.state, self.gas )
		elif scheme == 'AUSMDV':
			self.state = AUSMDV( self.domain, self.mesh, self.parameters, self.state, self.gas )
		t.toc('simulation time:')

		# post processing
		self.contourPanel.figure.clf()
		self.contourPanel.cax = self.contourPanel.figure.gca()
		plt.gca().xaxis.tick_bottom()
		self.contourPanel.cax.set_position([0.1, 0.16, 0.84, 0.82])
		#mpl.axes.Axes.clear(self.contourPanel.cax)
		cont = self.contourPanel.cax.contourf(self.mesh.xxc[1:-2,0:-1], self.mesh.yyc[1:-2,0:-1], \
							    		      self.state.Mach[1:-2,0:-1], 250, cmap=cm.jet)
		self.contourPanel.cax.axis('tight')
		self.contourPanel.cax.set_xlabel('x-coordinate (m)')
		self.contourPanel.cax.set_ylabel('y-coordinate (m)')
		self.contourPanel.canvas = FigureCanvas(self.contourPanel, -1, self.contourPanel.figure)

		# colorbar settings
		CB = self.contourPanel.figure.colorbar(cont, shrink=0.8, extend='both', ax=self.contourPanel.cax)
		CB.set_label('Mach Number', rotation=90)

		# residual plotting
		self.iterPanel.iax.plot(np.arange(1, len(self.state.res[0:self.state.n]), 1), self.state.res[1:self.state.n], linewidth=1)
		self.iterPanel.iax.set_xlabel('Iterations')
		self.iterPanel.iax.set_ylabel('Residual') 
		self.iterPanel.iax.get_lines()[0].set_color("black")
		self.iterPanel.iax.get_lines()[1].set_color("blue")
		self.iterPanel.iax.get_lines()[2].set_color("green")
		self.iterPanel.iax.get_lines()[3].set_color("red")
		self.iterPanel.iax.legend(['mdot', 'u', 'v', 'energy'], loc='center left', bbox_to_anchor=(1.05, 0.5))
		self.iterPanel.canvas = FigureCanvas(self.iterPanel, -1, self.iterPanel.figure)


		event.Skip()

	
	def call_grid( self, event ):
		import numpy as np
		from python.mesh.grid.gen_grid import mesh_wedge, mesh_airfoil
		from python.mesh.metrics.calc_cell_metrics import cellmetrics
		import matplotlib as mpl
		from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

		class domain:
			name = self.gridChoice.Strings[self.gridChoice.Selection]
			M = int(wx.grid.Grid.GetCellValue(self.domainGrid, 5, 0))
			N = int(wx.grid.Grid.GetCellValue(self.domainGrid, 6, 0))
			obj_start = float(wx.grid.Grid.GetCellValue(self.domainGrid, 2, 0))
			obj_end = float(wx.grid.Grid.GetCellValue(self.domainGrid, 3, 0))
			length = float(wx.grid.Grid.GetCellValue(self.domainGrid, 0, 0))
			height = float(wx.grid.Grid.GetCellValue(self.domainGrid, 1, 0))
			theta = np.deg2rad(float(wx.grid.Grid.GetCellValue(self.domainGrid, 4, 0)))

		if domain.name == "Wedge":
			xx, yy = mesh_wedge(domain)
		elif domain.name == "Airfoil":
			xx, yy = mesh_airfoil(domain)
		self.mesh = cellmetrics(xx, yy, domain)

		self.domain = domain

		# mesh plotting
		mpl.axes.Axes.clear(self.contourPanel.cax)
		self.contourPanel.cax.plot(self.mesh.xx, self.mesh.yy, color='blue', linewidth=0.5)
		self.contourPanel.cax.plot(np.transpose(self.mesh.xx), np.transpose(self.mesh.yy), color='blue', linewidth=0.5)
		self.contourPanel.cax.plot(self.mesh.xxc, self.mesh.yyc, 'gx', markersize=2)

		self.contourPanel.cax.set_xlim([np.min(self.mesh.xx[0,:]), np.max(self.mesh.xx[-1,:])])
		self.contourPanel.cax.set_ylim([np.min(self.mesh.yy[0,:]), domain.height])

		#self.contourPanel.cax.xaxis.tick_bottom()
		self.contourPanel.cax.set_xlabel('x-coordinate (m)')
		self.contourPanel.cax.set_ylabel('y-coordinate (m)')
		self.contourPanel.cax.axis('tight')
		self.contourPanel.canvas = FigureCanvas(self.contourPanel, -1, self.contourPanel.figure)

		sizer = wx.BoxSizer(wx.HORIZONTAL)
		sizer.Add(self.contourPanel.canvas, proportion=1, flag=wx.LEFT | wx.TOP | wx.GROW)
		self.SetSizer(sizer)

		event.Skip()

	
	def call_init( self, event ):
		import numpy as np
		from python.mesh.grid.gen_grid import mesh_wedge, mesh_airfoil
		from python.mesh.metrics.calc_cell_metrics import cellmetrics
		import matplotlib as mpl
		from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
		from matplotlib import cm
		from pytictoc import TicToc

		t = TicToc()

		# initialize state vector, simulation parameters and fluid properties
		class parameters:
			M_in = float(wx.grid.Grid.GetCellValue(self.parameterGrid, 0, 0))
			p_in = float(wx.grid.Grid.GetCellValue(self.parameterGrid, 1, 0))
			T_in = float(wx.grid.Grid.GetCellValue(self.parameterGrid, 2, 0))
			iterations = int(wx.grid.Grid.GetCellValue(self.simGrid, 1, 0))
			tolerance = float(wx.grid.Grid.GetCellValue(self.simGrid, 2, 0))
			CFL = float(wx.grid.Grid.GetCellValue(self.simGrid, 0, 0))
		class gas:
			gamma = 1.4
			Cp = 1006
			R = 287

		self.parameters = parameters
		self.gas = gas

		# initialize state vector, thermodynamic variables
		t.tic()
		from python.boundary.initialize import init_state
		self.state = init_state(self.domain, self.mesh, self.parameters, self.gas)
		t.toc('initialize time:')

		mpl.axes.Axes.clear(self.contourPanel.cax)
		cont = self.contourPanel.cax.contourf(self.mesh.xxc[1:-1,0:-1], self.mesh.yyc[1:-1,0:-1], \
							    		      self.state.V[1:-1,0:-1], 250, cmap=cm.jet)
		self.contourPanel.cax.axis('tight')
		# self.contourPanel.cax.set_xlim(np.min(self.mesh.xxc[1:-2,:]), np.max(self.mesh.xxc[1:-2,:]))
		# self.contourPanel.cax.set_ylim(np.min(self.mesh.yyc[:,1:-2]), np.max(self.mesh.yyc[:,1:-2]))
		self.contourPanel.cax.set_xlabel('x-coordinate (m)')
		self.contourPanel.cax.set_ylabel('y-coordinate (m)')
		self.contourPanel.canvas = FigureCanvas(self.contourPanel, -1, self.contourPanel.figure)

		event.Skip()
	

