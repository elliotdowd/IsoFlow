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

import subprocess
import sys

###########################################################################
## Class MainFrame
###########################################################################

class MainFrame ( wx.Frame ):
	
	def __init__( self, parent ):
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = wx.EmptyString, pos = wx.Point( 100,100 ), size = wx.Size( 740,800 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
		
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

		self.consolePanel = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.consolePanel.SetBackgroundColour( wx.Colour( 222, 222, 222 ) )
		MainSizer.Add( self.consolePanel, wx.GBPosition( 14, 0 ), wx.GBSpan( 5, 60 ), wx.ALL|wx.EXPAND, 5 )

		# TERMINAL COMMANDS
		self.consolePanel.command = wx.TextCtrl(self.consolePanel)
		self.consolePanel.result = wx.TextCtrl(self.consolePanel, style=wx.TE_MULTILINE)
		self.consolePanel.text = wx.TextCtrl(self.consolePanel, -1, style=wx.TE_MULTILINE|wx.TE_READONLY, size=(714, 60))

		redir=RedirectText(self.consolePanel.text)
		sys.stdout=redir

		
		self.iterPanel = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.iterPanel.SetBackgroundColour( wx.Colour( 222, 222, 222 ) )
		
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
		self.simGrid.SetRowLabelValue( 1, u"Iterations" )
		self.simGrid.SetRowLabelValue( 2, u"Residual Tol." )
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
		self.menuBar = wx.MenuBar( 0 )
		self.menuBar.SetFont( wx.Font( 12, 74, 90, 90, False, "Arial" ) )
		
		self.plotOptions = wx.Menu()
		self.contOptions = wx.Menu()
		self.mach = wx.MenuItem( self.contOptions, wx.ID_ANY, u"Mach Number", wx.EmptyString, wx.ITEM_NORMAL )
		self.contOptions.AppendItem( self.mach )
		
		self.pressure = wx.MenuItem( self.contOptions, wx.ID_ANY, u"Pressure", wx.EmptyString, wx.ITEM_NORMAL )
		self.contOptions.AppendItem( self.pressure )
		
		self.stagpressure = wx.MenuItem( self.contOptions, wx.ID_ANY, u"Stagnation Pressure", wx.EmptyString, wx.ITEM_NORMAL )
		self.contOptions.AppendItem( self.stagpressure )
		
		self.plotOptions.AppendSubMenu( self.contOptions, u"Contour" )
		
		self.cmOptions = wx.Menu()
		self.jet = wx.MenuItem( self.cmOptions, wx.ID_ANY, u"Jet", wx.EmptyString, wx.ITEM_NORMAL )
		self.cmOptions.AppendItem( self.jet )
		
		self.gray = wx.MenuItem( self.cmOptions, wx.ID_ANY, u"Greyscale", wx.EmptyString, wx.ITEM_NORMAL )
		self.cmOptions.AppendItem( self.gray )
		
		self.plotOptions.AppendSubMenu( self.cmOptions, u"Colormap" )
		
		self.menuBar.Append( self.plotOptions, u"Plot Options" ) 
		
		self.unitOptions = wx.Menu()
		self.metric2 = wx.MenuItem( self.unitOptions, wx.ID_ANY, u"Metric (kg-m-s-C)", wx.EmptyString, wx.ITEM_NORMAL )
		self.unitOptions.AppendItem( self.metric2 )
		
		self.metric1 = wx.MenuItem( self.unitOptions, wx.ID_ANY, u"Metric (kg-m-s-K)", wx.EmptyString, wx.ITEM_NORMAL )
		self.unitOptions.AppendItem( self.metric1 )
		
		self.menuBar.Append( self.unitOptions, u"Units" ) 
		
		self.SetMenuBar( self.menuBar )

		
		self.SetSizer( MainSizer )
		self.Layout()
		
		self.Centre( wx.BOTH )
		
		# Connect Events
		self.m_button3.Bind( wx.EVT_BUTTON, self.call_scheme )
		self.gridButton.Bind( wx.EVT_BUTTON, self.call_grid )
		self.initButton.Bind( wx.EVT_BUTTON, self.call_init )

		# initialize grid values
		self.init_grids()

	def __del__( self ):
		pass

	# initialize option grid values
	def init_grids( self ):
		# set domain row values
		self.domainGrid.SetCellValue( 0, 0, "1.5")
		self.domainGrid.SetCellValue( 1, 0, "1.3")
		self.domainGrid.SetCellValue( 2, 0, "0.5")
		self.domainGrid.SetCellValue( 3, 0, "1.1")
		self.domainGrid.SetCellValue( 4, 0, "20")
		self.domainGrid.SetCellValue( 5, 0, "48")
		self.domainGrid.SetCellValue( 6, 0, "36")

		# set initialization row values
		self.parameterGrid.SetCellValue( 0, 0, "3.0")
		self.parameterGrid.SetCellValue( 1, 0, "101325")
		self.parameterGrid.SetCellValue( 2, 0, "300")

		# set initialization row values
		self.simGrid.SetCellValue( 0, 0, "0.4")
		self.simGrid.SetCellValue( 1, 0, "1000")
		self.simGrid.SetCellValue( 2, 0, "-6")

	# Virtual event handlers, overide them in your derived class
	def call_grid( self, event ):
		import numpy as np
		from python.mesh.grid.gen_grid import mesh_wedge, mesh_airfoil
		from python.mesh.metrics.calc_cell_metrics import cellmetrics
		import matplotlib as mpl
		from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
		from pytictoc import TicToc

		t = TicToc()
		t.tic()

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

		print('________________________________________________________________________________________________________________________________________')
		print('Mesh elements: ' + str((domain.M+2) * (domain.N*2)))
		t.toc('Meshing time:')

		# mesh plotting
		self.contourPanel.figure.clf()
		self.contourPanel.cax = self.contourPanel.figure.gca()
		self.contourPanel.cax.set_position([0.1, 0.16, 0.84, 0.82])

		#mpl.axes.Axes.clear(self.contourPanel.cax)
		self.contourPanel.cax.plot(self.mesh.xx, self.mesh.yy, color='blue', linewidth=0.5)
		self.contourPanel.cax.plot(np.transpose(self.mesh.xx), np.transpose(self.mesh.yy), color='blue', linewidth=0.5)
		self.contourPanel.cax.plot(self.mesh.xxc, self.mesh.yyc, 'gx', markersize=2)

		self.contourPanel.cax.set_xlim([np.min(self.mesh.xx[0,:]), np.max(self.mesh.xx[-1,:])])
		self.contourPanel.cax.set_ylim([np.min(self.mesh.yy[0,:]), domain.height])

		#self.contourPanel.cax.xaxis.tick_bottom()
		self.contourPanel.cax.set_xlabel('x-coordinate (m)')
		self.contourPanel.cax.set_ylabel('y-coordinate (m)')
		self.contourPanel.cax.axis('equal')
		self.contourPanel.canvas = FigureCanvas(self.contourPanel, -1, self.contourPanel.figure)

		sizer = wx.BoxSizer(wx.HORIZONTAL)
		sizer.Add(self.contourPanel.canvas, proportion=1, flag=wx.LEFT | wx.TOP | wx.GROW)
		self.SetSizer(sizer)

		event.Skip()

	
	def call_init( self, event ):
		import numpy as np
		from python.mesh.grid.gen_grid import mesh_wedge, mesh_airfoil
		from python.mesh.metrics.calc_cell_metrics import cellmetrics
		from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
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

		print('________________________________________________________________________________________________________________________________________')
		t.toc('Initialize time:')

		self.call_contplot()

		event.Skip()
	

	def call_scheme( self, event ):

		from pytictoc import TicToc
		from python.finite_volume.AUSM.schemes import AUSM, AUSMplusup, AUSMDV

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

		self.call_contplot()
		self.call_resplot()

		event.Skip()


	def call_contplot(self):

		import numpy as np
		import matplotlib.pyplot as plt
		from matplotlib import cm
		import matplotlib as mpl
		from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

		# post processing
		self.contourPanel.figure.clf()
		self.contourPanel.cax = self.contourPanel.figure.gca()
		self.contourPanel.cax.set_position([0.12, 0.2, 0.84, 0.82])
		cont = self.contourPanel.cax.contourf(self.mesh.xxc[1:-2,1:-1], self.mesh.yyc[1:-2,1:-1], \
							    		      self.state.Mach[1:-2,1:-1], 250, cmap=cm.jet)
		self.contourPanel.cax.axis('equal')
		self.contourPanel.cax.set_xlabel('x-coordinate (m)')
		self.contourPanel.cax.set_ylabel('y-coordinate (m)')
		self.contourPanel.canvas = FigureCanvas(self.contourPanel, -1, self.contourPanel.figure)

		# colorbar settings
		ticks = np.linspace(round(np.min(self.state.Mach),2), round(np.max(self.state.Mach),2), 6)
		CB = self.contourPanel.figure.colorbar(cont, ticks=ticks, \
											   shrink=0.8, extend='both', ax=self.contourPanel.cax)
		CB.set_label('Mach Number', rotation=90)


	def call_resplot(self):

		import numpy as np
		import matplotlib.pyplot as plt
		from matplotlib import cm
		import matplotlib as mpl
		from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

		# residual plotting
		self.iterPanel.figure.clf()
		self.iterPanel.iax = self.iterPanel.figure.gca()
		self.iterPanel.iax.set_position([0.12, 0.06, 0.64, 0.72])
		self.iterPanel.iax.plot(np.arange(1, len(self.state.res[0:self.state.n]), 1), self.state.res[1:self.state.n], linewidth=1)
		self.iterPanel.iax.set_xlabel('Iterations')
		self.iterPanel.iax.set_ylabel('Residual') 
		self.iterPanel.iax.get_lines()[0].set_color("black")
		self.iterPanel.iax.get_lines()[1].set_color("blue")
		self.iterPanel.iax.get_lines()[2].set_color("green")
		self.iterPanel.iax.get_lines()[3].set_color("red")
		self.iterPanel.iax.legend([r"$\dot{m}$", 'u', 'v', r"$h_{t}$"], loc='center left', bbox_to_anchor=(1.025, 0.5), framealpha=0.0)
		self.iterPanel.canvas = FigureCanvas(self.iterPanel, -1, self.iterPanel.figure)

class RedirectText:
	def __init__(self,aWxTextCtrl):
		self.out=aWxTextCtrl

	def write(self,string):
		self.out.WriteText(string)
