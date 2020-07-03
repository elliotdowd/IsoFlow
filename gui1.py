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

import sys
from matplotlib import cm

###########################################################################
## Class MainFrame
###########################################################################

class MainFrame ( wx.Frame ):
	
	def __init__( self, parent ):
		wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = wx.EmptyString, pos = wx.Point( 100,100 ), size = wx.Size( 740,800 ), style=wx.DEFAULT_FRAME_STYLE ^ wx.RESIZE_BORDER )
		
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
		MainSizer.Add( self.contourPanel, wx.GBPosition( 1, 2 ), wx.GBSpan( 8, 54 ), wx.ALL|wx.EXPAND, 5 )
		
		self.consolePanel = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.consolePanel.SetBackgroundColour( wx.Colour( 222, 222, 222 ) )
		MainSizer.Add( self.consolePanel, wx.GBPosition( 14, 0 ), wx.GBSpan( 4, 60 ), wx.ALL|wx.EXPAND, 5 )

		# console redirecting
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

		self.gasOptions = wx.Menu()
		self.air = wx.MenuItem( self.gasOptions, 101, u"Air", wx.EmptyString, wx.ITEM_RADIO )
		self.gasOptions.Append( self.air )

		self.C02 = wx.MenuItem( self.gasOptions, 102, u"Carbon Dioxide", wx.EmptyString, wx.ITEM_RADIO )
		self.gasOptions.Append( self.C02 )

		self.H2 = wx.MenuItem( self.gasOptions, 103, u"Hydrogen", wx.EmptyString, wx.ITEM_RADIO )
		self.gasOptions.Append( self.H2 )
		
		self.gasOptions.AppendSeparator()
		
		self.thermalgas = wx.MenuItem( self.gasOptions, wx.ID_ANY, u"Thermally Perfect", wx.EmptyString, wx.ITEM_CHECK )
		self.gasOptions.Append( self.thermalgas )
		
		self.menuBar.Append( self.gasOptions, u"Gas" ) 
		
		self.plotOptions = wx.Menu()
		self.contOptions = wx.Menu()
		self.mach_change = wx.MenuItem( self.contOptions, 1, u"Mach Number", wx.EmptyString, wx.ITEM_RADIO )
		self.contOptions.Append( self.mach_change )

		self.velocity_change = wx.MenuItem( self.contOptions, 2, u"Velocity", wx.EmptyString, wx.ITEM_RADIO )
		self.contOptions.Append( self.velocity_change )

		self.quiver_change = wx.MenuItem( self.contOptions, 3, u"Velocity Streamlines", wx.EmptyString, wx.ITEM_RADIO )
		self.contOptions.Append( self.quiver_change )

		self.rho_change = wx.MenuItem( self.contOptions, 4, u"Density", wx.EmptyString, wx.ITEM_RADIO )
		self.contOptions.Append( self.rho_change )

		self.pressure_change = wx.MenuItem( self.contOptions, 5, u"Pressure", wx.EmptyString, wx.ITEM_RADIO )
		self.contOptions.Append( self.pressure_change )
		
		self.stagpressure_change = wx.MenuItem( self.contOptions, 6, u"Stagnation Pressure", wx.EmptyString, wx.ITEM_RADIO )
		self.contOptions.Append( self.stagpressure_change )
		
		self.temp_change = wx.MenuItem( self.contOptions, 7, u"Temperature", wx.EmptyString, wx.ITEM_RADIO )
		self.contOptions.Append( self.temp_change )
		
		self.stagtemp_change = wx.MenuItem( self.contOptions, 8, u"Stagnation Temperature", wx.EmptyString, wx.ITEM_RADIO )
		self.contOptions.Append( self.stagtemp_change )

		self.contOptions.AppendSeparator()

		self.coarse = wx.MenuItem( self.contOptions, wx.ID_ANY, u"Coarse", wx.EmptyString, wx.ITEM_RADIO )
		self.contOptions.Append( self.coarse )
		
		self.medium = wx.MenuItem( self.contOptions, wx.ID_ANY, u"Medium", wx.EmptyString, wx.ITEM_RADIO )
		self.contOptions.Append( self.medium )
		
		self.fine = wx.MenuItem( self.contOptions, wx.ID_ANY, u"Fine", wx.EmptyString, wx.ITEM_RADIO )
		self.contOptions.Append( self.fine )

		self.contOptions.AppendSeparator()
		
		self.label = wx.MenuItem( self.contOptions, wx.ID_ANY, u"Contour Labels", wx.EmptyString, wx.ITEM_CHECK )
		self.contOptions.Append( self.label )

		self.gradient = wx.MenuItem( self.contOptions, wx.ID_ANY, u"Gradient", wx.EmptyString, wx.ITEM_CHECK )
		self.contOptions.Append( self.gradient )
		
		self.plotOptions.AppendSubMenu( self.contOptions, u"Contour Options" )
		
		self.cmOptions = wx.Menu()
		self.jet = wx.MenuItem( self.cmOptions, wx.ID_ANY, u"Jet", wx.EmptyString, wx.ITEM_RADIO )
		self.cmOptions.Append( self.jet )

		self.magma = wx.MenuItem( self.cmOptions, wx.ID_ANY, u"Magma", wx.EmptyString, wx.ITEM_RADIO )
		self.cmOptions.Append( self.magma )
		
		self.gray = wx.MenuItem( self.cmOptions, wx.ID_ANY, u"Greyscale", wx.EmptyString, wx.ITEM_RADIO )
		self.cmOptions.Append( self.gray )
		
		self.plotOptions.AppendSubMenu( self.cmOptions, u"Colormap" )
		
		self.menuBar.Append( self.plotOptions, u"Plotting" ) 


		self.viewOptions = wx.Menu()
		self.expandCont = wx.MenuItem( self.viewOptions, wx.ID_ANY, u"Expand Contour Window", wx.EmptyString, wx.ITEM_NORMAL )
		self.viewOptions.Append( self.expandCont )

		self.axisOptions = wx.Menu()
		self.equal = wx.MenuItem( self.axisOptions, wx.ID_ANY, u"Equal", wx.EmptyString, wx.ITEM_RADIO )
		self.axisOptions.Append( self.equal )
		
		self.tight = wx.MenuItem( self.axisOptions, wx.ID_ANY, u"Tight", wx.EmptyString, wx.ITEM_RADIO )
		self.axisOptions.Append( self.tight )
		
		self.auto = wx.MenuItem( self.axisOptions, wx.ID_ANY, u"Auto", wx.EmptyString, wx.ITEM_RADIO )
		self.axisOptions.Append( self.auto )
		
		self.viewOptions.AppendSubMenu( self.axisOptions, u"Axis Options" )
		
		self.menuBar.Append( self.viewOptions, u"View" ) 

		
		self.unitOptions = wx.Menu()
		
		self.metric1 = wx.MenuItem( self.unitOptions, wx.ID_ANY, u"Metric (kg-m-s-K)", wx.EmptyString, wx.ITEM_RADIO )
		self.unitOptions.Append( self.metric1 )

		self.metric2 = wx.MenuItem( self.unitOptions, wx.ID_ANY, u"Metric (kg-m-s-°C)", wx.EmptyString, wx.ITEM_RADIO )
		self.unitOptions.Append( self.metric2 )

		self.imperial1 = wx.MenuItem( self.unitOptions, wx.ID_ANY, u"Imperial (lbm-ft-s-°F)", wx.EmptyString, wx.ITEM_RADIO )
		self.unitOptions.Append( self.imperial1 )

		self.imperial2 = wx.MenuItem( self.unitOptions, wx.ID_ANY, u"Imperial (slug-in-s-°R)", wx.EmptyString, wx.ITEM_RADIO )
		self.unitOptions.Append( self.imperial2 )
		
		self.menuBar.Append( self.unitOptions, u"Units" ) 
		
		self.SetMenuBar( self.menuBar )
		self.Centre( wx.BOTH )
		
		# Connect Events
		self.m_button3.Bind( wx.EVT_BUTTON, self.call_scheme )
		self.gridButton.Bind( wx.EVT_BUTTON, self.call_grid )
		self.initButton.Bind( wx.EVT_BUTTON, self.call_init )
		self.Bind( wx.EVT_MENU, self.gas_change, id = self.air.GetId() )
		self.Bind( wx.EVT_MENU, self.thermalgas_change, id = self.thermalgas.GetId() )
		self.Bind( wx.EVT_MENU, self.mach, id = self.mach_change.GetId() )
		self.Bind( wx.EVT_MENU, self.velocity, id = self.velocity_change.GetId() )
		self.Bind( wx.EVT_MENU, self.quiver, id = self.quiver_change.GetId() )
		self.Bind( wx.EVT_MENU, self.rho, id = self.rho_change.GetId() )
		self.Bind( wx.EVT_MENU, self.pressure, id = self.pressure_change.GetId() )
		self.Bind( wx.EVT_MENU, self.stagp, id = self.stagpressure_change.GetId() )
		self.Bind( wx.EVT_MENU, self.temp, id = self.temp_change.GetId() )
		self.Bind( wx.EVT_MENU, self.stagtemp, id = self.stagtemp_change.GetId() )
		self.Bind( wx.EVT_MENU, self.coarse_change, id = self.coarse.GetId() )
		self.Bind( wx.EVT_MENU, self.medium_change, id = self.medium.GetId() )
		self.Bind( wx.EVT_MENU, self.fine_change, id = self.fine.GetId() )
		self.Bind( wx.EVT_MENU, self.label_change, id = self.label.GetId() )
		self.Bind( wx.EVT_MENU, self.gradient_change, id = self.gradient.GetId() )
		self.Bind( wx.EVT_MENU, self.jet_change, id = self.jet.GetId() )
		self.Bind( wx.EVT_MENU, self.magma_change, id = self.magma.GetId() )
		self.Bind( wx.EVT_MENU, self.gray_change, id = self.gray.GetId() )
		self.Bind( wx.EVT_MENU, self.metric1_change, id = self.metric1.GetId())
		self.Bind( wx.EVT_MENU, self.metric2_change, id = self.metric2.GetId())
		self.Bind( wx.EVT_MENU, self.imperial1_change, id = self.imperial1.GetId())
		self.Bind( wx.EVT_MENU, self.imperial2_change, id = self.imperial2.GetId())
		self.Bind( wx.EVT_MENU, self.expandWindow, id = self.expandCont.GetId() )
		self.Bind( wx.EVT_MENU, self.equal_change, id = self.equal.GetId() )
		self.Bind( wx.EVT_MENU, self.tight_change, id = self.tight.GetId() )
		self.Bind( wx.EVT_MENU, self.auto_change, id = self.auto.GetId() )


		# initialize grid values and class attributes
		self.init_grids()
		self.contQuantity = 'Mach'
		self.cmOption = cm.jet
		class units:
			mass = 'kg'
			def conv_mass(m):
				conv = m
				return conv
			length = 'm'
			def conv_length(l):
				conv = l
				return conv
			time = 's'
			def conv_time(t):
				conv = t
				return conv
			temp = 'K'
			def conv_temp(T):
				conv = T
				return conv
			press = 'kPa'
			def conv_press(p):
				conv = p / 1000
				return conv
		self.units = units
		self.axisOption = 'equal'
		self.contGrad = 8
		self.labeled = False
		self.gradient = ''
		self.gasSelect = 'Air'
		self.thermoModel = 'cpg'



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

		# length unit conversion
		cl = self.units.conv_length(1)

		class domain:
			name = self.gridChoice.Strings[self.gridChoice.Selection]
			M = int(wx.grid.Grid.GetCellValue(self.domainGrid, 5, 0))
			N = int(wx.grid.Grid.GetCellValue(self.domainGrid, 6, 0))
			obj_start = float(wx.grid.Grid.GetCellValue(self.domainGrid, 2, 0)) / cl
			obj_end = float(wx.grid.Grid.GetCellValue(self.domainGrid, 3, 0)) / cl
			length = float(wx.grid.Grid.GetCellValue(self.domainGrid, 0, 0)) / cl
			height = float(wx.grid.Grid.GetCellValue(self.domainGrid, 1, 0)) / cl
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
		self.contourPanel.cax.plot(self.mesh.xx * cl, self.mesh.yy * cl, color='blue', linewidth=0.5)
		self.contourPanel.cax.plot(np.transpose(self.mesh.xx) * cl, np.transpose(self.mesh.yy) * cl, color='blue', linewidth=0.5)
		self.contourPanel.cax.plot(self.mesh.xxc * cl, self.mesh.yyc * cl, 'gx', markersize=2)

		self.contourPanel.cax.set_xlim([np.min(self.mesh.xx[0,:]), np.max(self.mesh.xx[-1,:])])
		self.contourPanel.cax.set_ylim([np.min(self.mesh.yy[0,:]), domain.height])

		#self.contourPanel.cax.xaxis.tick_bottom()
		self.contourPanel.cax.set_xlabel('x-coordinate ' + '(' + self.units.length + ')')
		self.contourPanel.cax.set_ylabel('y-coordinate ' + '(' + self.units.length + ')')
		self.contourPanel.cax.axis(self.axisOption)
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
		import python.finite_volume.gasdata as gasdata

		t = TicToc()

		# initialize state vector, simulation parameters and fluid properties
		class parameters:
			M_in = float(wx.grid.Grid.GetCellValue(self.parameterGrid, 0, 0))
			p_in = float(wx.grid.Grid.GetCellValue(self.parameterGrid, 1, 0))
			T_in = float(wx.grid.Grid.GetCellValue(self.parameterGrid, 2, 0))
			iterations = int(wx.grid.Grid.GetCellValue(self.simGrid, 1, 0))
			tolerance = float(wx.grid.Grid.GetCellValue(self.simGrid, 2, 0))
			CFL = float(wx.grid.Grid.GetCellValue(self.simGrid, 0, 0))

		self.parameters = parameters
		self.gas = gasdata.air_tpg

		# initialize state vector, thermodynamic variables
		t.tic()
		from python.boundary.initialize import init_state
		self.state = init_state(self.domain, self.mesh, self.parameters, self.gas)

		print('________________________________________________________________________________________________________________________________________')
		t.toc('Initialize time:')

		self.call_contplot(self.contourPanel)

		event.Skip()

	def call_scheme( self, event ):

		from pytictoc import TicToc
		from python.finite_volume.AUSM.schemes import AUSM, AUSMplusup, AUSMDV
		import python.finite_volume.gasdata as gasdata

		t = TicToc()

		# run AUSM family scheme
		t.tic()
		scheme = self.schemeChoice.Strings[self.schemeChoice.Selection]

		class parameters:
			M_in = float(wx.grid.Grid.GetCellValue(self.parameterGrid, 0, 0))
			p_in = float(wx.grid.Grid.GetCellValue(self.parameterGrid, 1, 0))
			T_in = float(wx.grid.Grid.GetCellValue(self.parameterGrid, 2, 0))
			iterations = int(wx.grid.Grid.GetCellValue(self.simGrid, 1, 0))
			tolerance = float(wx.grid.Grid.GetCellValue(self.simGrid, 2, 0))
			CFL = float(wx.grid.Grid.GetCellValue(self.simGrid, 0, 0))

		self.parameters = parameters

		if self.thermoModel == 'cpg':
			if self.air.IsChecked() == True:
				self.gas = gasdata.air_cpg
			elif self.C02.IsChecked() == True:
				self.gas = gasdata.C02_cpg
			elif self.H2.IsChecked() == True:
				self.gas = gasdata.H2_cpg
		elif self.thermoModel == 'tpg':
			if self.air.IsChecked() == True:
				self.gas = gasdata.air_tpg
			elif self.C02.IsChecked() == True:
				self.gas = gasdata.C02_tpg
			elif self.H2.IsChecked() == True:
				self.gas = gasdata.H2_tpg
		
		if scheme == 'AUSM':
			self.state = AUSM( self.domain, self.mesh, self.parameters, self.state, self.gas )
		elif scheme == 'AUSM+up':
			self.state = AUSMplusup( self.domain, self.mesh, self.parameters, self.state, self.gas )
		elif scheme == 'AUSMDV':
			self.state = AUSMDV( self.domain, self.mesh, self.parameters, self.state, self.gas )
		t.toc('simulation time:')

		self.call_contplot(self.contourPanel)
		self.call_resplot()

		event.Skip()

	def call_contplot(self, panel):

		import numpy as np
		import matplotlib.pyplot as plt
		from matplotlib import cm
		import matplotlib as mpl
		from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

		# panel input as self.contourPanel

		# post processing
		panel.figure.clf()
		panel.cax = panel.figure.gca()
		panel.cax.set_facecolor((0.4, 0.4, 0.4))
		panel.cax.set_position([0.12, 0.2, 0.84, 0.82])

		contQuantity = self.contQuantity + ' ' + self.gradient
		cl = self.units.conv_length(1)

		if contQuantity == 'Mach ':
			cont = panel.cax.contourf(cl*self.mesh.xxc[0:-1,1:-1], cl*self.mesh.yyc[0:-1,1:-1], \
							    		      	  self.state.Mach[0:-1,1:-1], self.contGrad, cmap=self.cmOption)
			# colorbar settings
			ticks = np.linspace(round(np.min(self.state.Mach),2), round(np.max(self.state.Mach),2), 6)
			CB = panel.figure.colorbar(cont, ticks=ticks, \
												shrink=0.8, extend='both', ax=panel.cax)
			CB.set_label(contQuantity, rotation=90)
		elif contQuantity == 'Velocity ':
			velocity = (cl/self.units.conv_time(1)) * self.state.vel[0:-1,1:-1]
			cont = panel.cax.contourf(cl*self.mesh.xxc[0:-1,1:-1], cl*self.mesh.yyc[0:-1,1:-1], \
							    		      	  velocity, self.contGrad, cmap=self.cmOption)
			# colorbar settings
			ticks = np.linspace(round(np.min(velocity),0), round(np.max(velocity),0), 6)
			CB = panel.figure.colorbar(cont, ticks=ticks, \
												shrink=0.8, extend='both', ax=panel.cax)
			CB.set_label(contQuantity + ' (' + self.units.length + '/' + self.units.time + ')', rotation=90)
		elif contQuantity == 'Velocity Quiver ':
			cont = panel.cax.quiver(cl*self.mesh.xxc[0:-1,1:-1], cl*self.mesh.yyc[0:-1,1:-1], \
								  				self.state.u[0:-1,1:-1], self.state.v[0:-1,1:-1], \
												self.state.vel[0:-1,1:-1], cmap=self.cmOption, pivot='tip', \
												angles='uv', scale_units='width', scale=0.8*self.domain.M*np.max(self.state.vel[0:-1,1:-1]))
			# colorbar settings
			ticks = np.linspace(round(np.min(self.state.vel),0), round(np.max(self.state.vel),0), 6)
			CB = panel.figure.colorbar(cont, ticks=ticks, \
												shrink=0.8, extend='both', ax=panel.cax)
			CB.set_label(contQuantity + ' (' + self.units.length + '/' + self.units.time + ')', rotation=90)
		elif contQuantity == 'Density ':
			rho = self.units.conv_mass(1)/self.units.conv_length(1)**3
			cont = panel.cax.contourf(cl*self.mesh.xxc[0:-1,1:-1], cl*self.mesh.yyc[0:-1,1:-1], \
							    		      	  rho*self.state.Q[0:-1,1:-1,0], self.contGrad, cmap=self.cmOption)
			# colorbar settings
			ticks = np.linspace(round(np.min(rho*self.state.Q[0:-1,1:-1,0]),9), round(np.max(rho*self.state.Q[0:-1,1:-1,0]),9), 6)
			CB = panel.figure.colorbar(cont, ticks=ticks, \
												shrink=0.8, extend='both', ax=panel.cax)
			CB.set_label(contQuantity + ' (' + self.units.mass + '/' + self.units.length + '$^3$' + ')', rotation=90)
		elif contQuantity == 'Density Gradient':
			from python.finite_volume.helper import grad
			rho = self.units.conv_mass(1)/self.units.conv_length(1)**4
			rhograd = grad(cl*self.mesh.xxc, cl*self.mesh.yyc, rho*self.state.Q[:,:,0])
			cont = panel.cax.contourf(cl*self.mesh.xxc[1:-1,1:-1], cl*self.mesh.yyc[1:-1,1:-1], \
							    		      	  rhograd, self.contGrad, cmap=self.cmOption)
			# colorbar settings
			ticks = np.linspace(round(np.min(rhograd), 9), round(np.max(rhograd), 9), 6)
			CB = panel.figure.colorbar(cont, ticks=ticks, \
												shrink=0.8, extend='both', ax=panel.cax)
			CB.set_label(contQuantity + ' (' + self.units.mass + '/' + self.units.length + '$^4$' + ')', rotation=90)
		elif contQuantity == 'Pressure ':
			pressure = self.units.conv_press(self.state.p[0:-1,1:-1])
			cont = panel.cax.contourf(cl*self.mesh.xxc[0:-1,1:-1], cl*self.mesh.yyc[0:-1,1:-1], \
							    		      	  pressure, self.contGrad, cmap=self.cmOption)
			# colorbar settings
			ticks = np.linspace(round(np.min(pressure),0), round(np.max(pressure),0), 6)
			CB = panel.figure.colorbar(cont, ticks=ticks, \
												shrink=0.8, extend='both', ax=panel.cax)
			CB.set_label(contQuantity + ' (' + self.units.press + ')', rotation=90)
		elif contQuantity == 'Pressure Gradient':
			from python.finite_volume.helper import grad
			pg = self.units.conv_press(1) / self.units.conv_length(1)
			pgrad = grad(cl*self.mesh.xxc, cl*self.mesh.yyc, pg*self.state.p)
			cont = panel.cax.contourf(cl*self.mesh.xxc[1:-1,1:-1], cl*self.mesh.yyc[1:-1,1:-1], \
							    		      	  pgrad, self.contGrad, cmap=self.cmOption)
			# colorbar settings
			ticks = np.linspace(round(np.min(pgrad), 9), round(np.max(pgrad), 9), 6)
			CB = panel.figure.colorbar(cont, ticks=ticks, \
												shrink=0.8, extend='both', ax=panel.cax)
			CB.set_label(contQuantity + ' (' + self.units.press + '/' + self.units.length + ')', rotation=90)
		elif contQuantity == 'Stagnation Pressure ':
			p0 = self.units.conv_press(self.state.p0[0:-1,1:-1])
			cont = self.contourPanel.cax.contourf(cl*self.mesh.xxc[0:-1,1:-1], cl*self.mesh.yyc[0:-1,1:-1], \
							    		      	  p0, self.contGrad, cmap=self.cmOption)
			# colorbar settings
			ticks = np.linspace(round(np.min(p0),0), round(np.max(p0),0), 6)
			CB = self.contourPanel.figure.colorbar(cont, ticks=ticks, \
												shrink=0.8, extend='both', ax=panel.cax)
			CB.set_label(contQuantity + ' (' + self.units.press + ')', rotation=90)
		elif contQuantity == 'Temperature ':
			temperature = self.units.conv_temp(self.state.T[0:-1,1:-1])
			cont = panel.cax.contourf(cl*self.mesh.xxc[0:-1,1:-1], cl*self.mesh.yyc[0:-1,1:-1], \
							    		      	  temperature, self.contGrad, cmap=self.cmOption)
			# colorbar settings
			ticks = np.linspace(round(np.min(temperature),0), round(np.max(temperature),0), 6)
			CB = panel.figure.colorbar(cont, ticks=ticks, \
												shrink=0.8, extend='both', ax=panel.cax)
			CB.set_label(contQuantity + ' (' + self.units.temp + ')', rotation=90)
		elif contQuantity == 'Temperature Gradient':
			from python.finite_volume.helper import grad
			ct = self.units.conv_temp(1) / self.units.conv_length(1)
			tgrad = grad(cl*self.mesh.xxc, cl*self.mesh.yyc, ct*self.state.T)
			cont = panel.cax.contourf(cl*self.mesh.xxc[1:-1,1:-1], cl*self.mesh.yyc[1:-1,1:-1], \
							    		      	  tgrad, self.contGrad, cmap=self.cmOption)
			# colorbar settings
			ticks = np.linspace(round(np.min(tgrad), 0), round(np.max(tgrad), 0), 6)
			CB = panel.figure.colorbar(cont, ticks=ticks, \
												shrink=0.8, extend='both', ax=panel.cax)
			CB.set_label(contQuantity + ' (' + self.units.temp + '/' + self.units.length + ')', rotation=90)
		elif contQuantity == 'Stagnation Temperature ':
			T0 = self.units.conv_temp(self.state.T0[0:-1,1:-1])
			cont = panel.cax.contourf(cl*self.mesh.xxc[0:-1,1:-1], cl*self.mesh.yyc[0:-1,1:-1], \
							    		      	  T0, self.contGrad, cmap=self.cmOption)
			# colorbar settings
			ticks = np.linspace(round(np.min(T0),0), round(np.max(T0),0), 6)
			CB = panel.figure.colorbar(cont, ticks=ticks, \
												shrink=0.8, extend='both', ax=panel.cax)
			CB.set_label(contQuantity + ' (' + self.units.temp + ')', rotation=90)

		# set up contour labels
		if self.contQuantity != 'Velocity Quiver':
			if self.labeled:
				if len(cont.levels) > self.contGrad**(1/3):
					self.contourPanel.cax.clabel(cont, cont.levels[0:self.contGrad:int(self.contGrad*0.5/self.contGrad**(1/3))], fmt='%2.3f', colors='w', fontsize=8)
				else:
					self.contourPanel.cax.clabel(cont, fmt='%2.3f', colors='w', fontsize=8)

		# plot settings
		panel.cax.xaxis.tick_bottom()
		panel.cax.set_xlabel('x-coordinate' + ' (' + self.units.length + ')')
		panel.cax.set_ylabel('y-coordinate' + ' (' + self.units.length + ')')
		panel.cax.axis(self.axisOption)
		panel.canvas = FigureCanvas(panel, -1, panel.figure)

	def call_resplot(self):

		import numpy as np
		import matplotlib.pyplot as plt
		from matplotlib import cm
		import matplotlib as mpl
		from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

		# residual plotting
		self.iterPanel.figure.clf()
		self.iterPanel.iax = self.iterPanel.figure.gca()
		self.iterPanel.iax.set_position([0.14, 0.06, 0.68, 0.72])
		self.iterPanel.iax.plot(np.arange(1, len(self.state.res[0:self.state.n]), 1), self.state.res[1:self.state.n], linewidth=0.8)
		self.iterPanel.iax.set_xlabel('Iterations')
		self.iterPanel.iax.set_ylabel('Residual') 
		self.iterPanel.iax.get_lines()[0].set_color("black")
		self.iterPanel.iax.get_lines()[1].set_color("blue")
		self.iterPanel.iax.get_lines()[2].set_color("green")
		self.iterPanel.iax.get_lines()[3].set_color("red")
		self.iterPanel.iax.legend([r"$\dot{m}$", 'u', 'v', r"$h_{t}$"], loc='center left', bbox_to_anchor=(1.025, 0.5), framealpha=0.0)
		self.iterPanel.canvas = FigureCanvas(self.iterPanel, -1, self.iterPanel.figure)


	# menubar events
	def mach( self, event ):
		self.contQuantity = 'Mach'
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()

	def velocity( self, event ):
		self.contQuantity = 'Velocity'
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()

	def quiver( self, event ):
		self.contQuantity = 'Velocity Quiver'
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()
	
	def rho( self, event ):
		self.contQuantity = 'Density'
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()

	def pressure( self, event ):
		self.contQuantity = 'Pressure'
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()
	
	def stagp( self, event ):
		self.contQuantity = 'Stagnation Pressure'
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()

	def temp( self, event ):
		self.contQuantity = 'Temperature'
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()
	
	def stagtemp( self, event ):
		self.contQuantity = 'Stagnation Temperature'
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()

	def jet_change( self, event ):
		from matplotlib import cm
		self.cmOption = cm.jet
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()

	def magma_change( self, event ):
		from matplotlib import cm
		self.cmOption = cm.magma
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()

	def gray_change( self, event ):
		from matplotlib import cm
		self.cmOption = cm.gray
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()

	def metric1_change( self, event ):
		class units:
			mass = 'kg'
			def conv_mass(m):
				conv = m
				return conv
			length = 'm'
			def conv_length(l):
				conv = l
				return conv
			time = 's'
			def conv_time(t):
				conv = t
				return conv
			temp = 'K'
			def conv_temp(T):
				conv = T
				return conv
			press = 'kPa'
			def conv_press(p):
				conv = p / 1000
				return conv
		self.units = units
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)

	def metric2_change( self, event ):
		class units:
			mass = 'kg'
			def conv_mass(m):
				conv = m
				return conv
			length = 'm'
			def conv_length(l):
				conv = l
				return conv
			time = 's'
			def conv_time(t):
				conv = t
				return conv
			temp = '°C'
			def conv_temp(T):
				conv = T - 272.15
				return conv
			press = 'kPa'
			def conv_press(p):
				conv = p / 1000
				return conv
		self.units = units
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)

	def imperial1_change( self, event ):
		class units:
			mass = 'lbm'
			def conv_mass(m):
				conv = m * 2.20462
				return conv
			length = 'ft'
			def conv_length(l):
				conv = l * 3.28084
				return conv
			time = 's'
			def conv_time(t):
				conv = t
				return conv
			temp = '°F'
			def conv_temp(T):
				conv = (T - 272.15) * 9/5 + 32
				return conv
			press = 'psi'
			def conv_press(p):
				conv = p * 0.000145038
				return conv
		self.units = units
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)

	def imperial2_change( self, event ):
		class units:
			mass = 'slug'
			def conv_mass(m):
				conv = m * 0.0685218
				return conv
			length = 'in'
			def conv_length(l):
				conv = l * 39.3701
				return conv
			time = 's'
			def conv_time(t):
				conv = t
				return conv
			temp = '°R'
			def conv_temp(T):
				conv = (T - 272.15) * 9/5 + 32 + 491.67
				return conv
			press = 'lbf/ft$^2$'
			def conv_press(p):
				conv = p * 0.02088545226628
				return conv
		self.units = units
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)

	def equal_change( self, event ):
		self.axisOption = 'equal'
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()
	
	def tight_change( self, event ):
		self.axisOption = 'tight'
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()
	
	def auto_change( self, event ):
		self.axisOption = 'auto'
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()

	def coarse_change( self, event ):
		self.contGrad = 8
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()
	
	def medium_change( self, event ):
		self.contGrad = 64
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()
	
	def fine_change( self, event ):
		self.contGrad = 512
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()

	def label_change( self, event ):
		if self.labeled == False:
			self.labeled = True
		else:
			self.labeled = False
		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()

	def gradient_change( self, event ):
		if self.gradient == 'Gradient':
			self.gradient = ''
			wx.MenuBar.Enable(self.menuBar, 1, True)
			wx.MenuBar.Enable(self.menuBar, 2, True)
			wx.MenuBar.Enable(self.menuBar, 3, True)

			wx.MenuBar.Enable(self.menuBar, 6, True)
			#wx.MenuBar.Enable(self.menuBar, 7, True)
			wx.MenuBar.Enable(self.menuBar, 8, True)

		else:
			self.gradient = 'Gradient'
			wx.MenuBar.Enable(self.menuBar, 1, False)
			wx.MenuBar.Enable(self.menuBar, 2, False)
			wx.MenuBar.Enable(self.menuBar, 3, False)

			wx.MenuBar.Enable(self.menuBar, 6, False)
			#wx.MenuBar.Enable(self.menuBar, 7, False)
			wx.MenuBar.Enable(self.menuBar, 8, False)

		if hasattr(self, 'state'):
			self.call_contplot(self.contourPanel)
		event.Skip()

	def gas_change( self, event ):
		event.Skip()
	
	def thermalgas_change( self, event ):
		if self.thermoModel == 'tpg':
			self.thermoModel = 'cpg'
		else: 
			self.thermoModel = 'tpg'
		event.Skip()

	# open contour plot in new window
	def expandWindow( self, event ):
		#import matplotlib.pyplot as plt
		self.new = NewWindow(parent=None)
		self.call_contplot(self.new.contourPanel)
		self.new.Show()
		event.Skip()

	# open informational window for gases

class RedirectText:
	def __init__(self,aWxTextCtrl):
		self.out=aWxTextCtrl

	def write(self,string):
		self.out.WriteText(string)

		# redir=RedirectText(self.consolePanel.text)
		# sys.stdout=redir


class NewWindow(wx.Frame):
	def __init__(self, parent):
		import matplotlib.pyplot as plt
		import numpy as np
		wx.Frame.__init__( self, parent, title = 'Fullscreen Contour Plot',\
						   size = wx.Size( 1020,720 ), style=wx.DEFAULT_FRAME_STYLE )

		self.SetBackgroundColour( wx.Colour( 256, 256, 256 ) )

		self.contourPanel = wx.Panel( self, wx.ID_ANY, wx.DefaultPosition, wx.DefaultSize, wx.TAB_TRAVERSAL )
		self.contourPanel.SetBackgroundColour( wx.Colour( 256, 256, 256 ) )

		self.contourPanel.figure = plt.figure( dpi=100, figsize=(10, 10/1.4473))
		self.contourPanel.cax = self.contourPanel.figure.gca()
		self.contourPanel.cax.set_facecolor((0.4, 0.4, 0.4))
		self.contourPanel.cax.set_position([0.15, 0.1, 0.8, 0.76])

