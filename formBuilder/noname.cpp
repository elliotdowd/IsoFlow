///////////////////////////////////////////////////////////////////////////
// C++ code generated with wxFormBuilder (version Jun 17 2015)
// http://www.wxformbuilder.org/
//
// PLEASE DO "NOT" EDIT THIS FILE!
///////////////////////////////////////////////////////////////////////////

#include "noname.h"

///////////////////////////////////////////////////////////////////////////

MainFrame::MainFrame( wxWindow* parent, wxWindowID id, const wxString& title, const wxPoint& pos, const wxSize& size, long style ) : wxFrame( parent, id, title, pos, size, style )
{
	this->SetSizeHints( wxDefaultSize, wxDefaultSize );
	
	wxGridBagSizer* MainSizer;
	MainSizer = new wxGridBagSizer( 0, 0 );
	MainSizer->SetFlexibleDirection( wxBOTH );
	MainSizer->SetNonFlexibleGrowMode( wxFLEX_GROWMODE_SPECIFIED );
	
	domainGrid = new wxGrid( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, 0 );
	
	// Grid
	domainGrid->CreateGrid( 6, 1 );
	domainGrid->EnableEditing( true );
	domainGrid->EnableGridLines( true );
	domainGrid->EnableDragGridSize( false );
	domainGrid->SetMargins( 0, 0 );
	
	// Columns
	domainGrid->SetColSize( 0, 72 );
	domainGrid->EnableDragColMove( false );
	domainGrid->EnableDragColSize( true );
	domainGrid->SetColLabelSize( 0 );
	domainGrid->SetColLabelAlignment( wxALIGN_CENTRE, wxALIGN_CENTRE );
	
	// Rows
	domainGrid->EnableDragRowSize( true );
	domainGrid->SetRowLabelSize( 102 );
	domainGrid->SetRowLabelValue( 0, wxT("Length") );
	domainGrid->SetRowLabelValue( 1, wxT("Height") );
	domainGrid->SetRowLabelValue( 2, wxT("Wedge Start") );
	domainGrid->SetRowLabelValue( 3, wxT("Wedge Angle") );
	domainGrid->SetRowLabelValue( 4, wxT("Horiz. Cells") );
	domainGrid->SetRowLabelValue( 5, wxT("Vertical Cells") );
	domainGrid->SetRowLabelAlignment( wxALIGN_CENTRE, wxALIGN_CENTRE );
	
	// Label Appearance
	
	// Cell Defaults
	domainGrid->SetDefaultCellAlignment( wxALIGN_LEFT, wxALIGN_TOP );
	MainSizer->Add( domainGrid, wxGBPosition( 1, 0 ), wxGBSpan( 1, 1 ), wxALL|wxEXPAND, 5 );
	
	contourPanel = new wxPanel( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxTAB_TRAVERSAL );
	MainSizer->Add( contourPanel, wxGBPosition( 1, 2 ), wxGBSpan( 7, 54 ), wxALL|wxEXPAND, 5 );
	
	parameterGrid = new wxGrid( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, 0 );
	
	// Grid
	parameterGrid->CreateGrid( 3, 1 );
	parameterGrid->EnableEditing( true );
	parameterGrid->EnableGridLines( true );
	parameterGrid->EnableDragGridSize( false );
	parameterGrid->SetMargins( 0, 0 );
	
	// Columns
	parameterGrid->SetColSize( 0, 54 );
	parameterGrid->EnableDragColMove( false );
	parameterGrid->EnableDragColSize( true );
	parameterGrid->SetColLabelSize( 0 );
	parameterGrid->SetColLabelAlignment( wxALIGN_CENTRE, wxALIGN_CENTRE );
	
	// Rows
	parameterGrid->EnableDragRowSize( true );
	parameterGrid->SetRowLabelSize( 120 );
	parameterGrid->SetRowLabelValue( 0, wxT("Inlet Mach #") );
	parameterGrid->SetRowLabelValue( 1, wxT("Inlet Pres. (Pa)") );
	parameterGrid->SetRowLabelValue( 2, wxT("Inlet Temp. (K)") );
	parameterGrid->SetRowLabelAlignment( wxALIGN_CENTRE, wxALIGN_CENTRE );
	
	// Label Appearance
	
	// Cell Defaults
	parameterGrid->SetDefaultCellAlignment( wxALIGN_LEFT, wxALIGN_TOP );
	MainSizer->Add( parameterGrid, wxGBPosition( 6, 0 ), wxGBSpan( 1, 1 ), wxALL|wxEXPAND, 5 );
	
	simGrid = new wxGrid( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, 0 );
	
	// Grid
	simGrid->CreateGrid( 3, 1 );
	simGrid->EnableEditing( true );
	simGrid->EnableGridLines( true );
	simGrid->EnableDragGridSize( false );
	simGrid->SetMargins( 0, 0 );
	
	// Columns
	simGrid->SetColSize( 0, 54 );
	simGrid->EnableDragColMove( false );
	simGrid->EnableDragColSize( true );
	simGrid->SetColLabelSize( 0 );
	simGrid->SetColLabelAlignment( wxALIGN_CENTRE, wxALIGN_CENTRE );
	
	// Rows
	simGrid->EnableDragRowSize( true );
	simGrid->SetRowLabelSize( 120 );
	simGrid->SetRowLabelValue( 0, wxT("Max. CFL") );
	simGrid->SetRowLabelValue( 1, wxT("Inlet Pres. (Pa)") );
	simGrid->SetRowLabelValue( 2, wxT("Residual Tol.") );
	simGrid->SetRowLabelValue( 3, wxT("Iterations") );
	simGrid->SetRowLabelValue( 4, wxEmptyString );
	simGrid->SetRowLabelAlignment( wxALIGN_CENTRE, wxALIGN_CENTRE );
	
	// Label Appearance
	
	// Cell Defaults
	simGrid->SetDefaultCellAlignment( wxALIGN_LEFT, wxALIGN_TOP );
	MainSizer->Add( simGrid, wxGBPosition( 10, 0 ), wxGBSpan( 1, 1 ), wxALL, 5 );
	
	m_staticText11 = new wxStaticText( this, wxID_ANY, wxT("Inlet Parameters"), wxDefaultPosition, wxDefaultSize, wxALIGN_CENTRE );
	m_staticText11->Wrap( -1 );
	m_staticText11->SetFont( wxFont( 10, 74, 90, 92, false, wxT("Arial") ) );
	
	MainSizer->Add( m_staticText11, wxGBPosition( 5, 0 ), wxGBSpan( 1, 1 ), wxALL|wxALIGN_CENTER_HORIZONTAL, 5 );
	
	m_staticText111 = new wxStaticText( this, wxID_ANY, wxT("Simulation Options"), wxDefaultPosition, wxDefaultSize, wxALIGN_CENTRE );
	m_staticText111->Wrap( -1 );
	m_staticText111->SetFont( wxFont( 10, 74, 90, 92, false, wxT("Arial") ) );
	
	MainSizer->Add( m_staticText111, wxGBPosition( 9, 0 ), wxGBSpan( 1, 1 ), wxALL|wxALIGN_CENTER_HORIZONTAL, 5 );
	
	meshType = new wxListBox( this, wxID_ANY, wxDefaultPosition, wxDefaultSize, 0, NULL, wxLB_SINGLE );
	meshType->Append( wxT("Wedge") );
	meshType->Append( wxT("Airfoil") );
	MainSizer->Add( meshType, wxGBPosition( 2, 0 ), wxGBSpan( 1, 1 ), wxALL|wxALIGN_RIGHT|wxEXPAND, 5 );
	
	m_button3 = new wxButton( this, wxID_ANY, wxT("Run Simulation"), wxDefaultPosition, wxDefaultSize, 0 );
	MainSizer->Add( m_button3, wxGBPosition( 11, 0 ), wxGBSpan( 1, 1 ), wxALL|wxEXPAND, 5 );
	
	gridButton = new wxButton( this, wxID_ANY, wxT("Generate Grid"), wxDefaultPosition, wxDefaultSize, 0 );
	MainSizer->Add( gridButton, wxGBPosition( 3, 0 ), wxGBSpan( 1, 1 ), wxALL|wxEXPAND, 5 );
	
	initButton = new wxButton( this, wxID_ANY, wxT("Initialize"), wxDefaultPosition, wxDefaultSize, 0 );
	MainSizer->Add( initButton, wxGBPosition( 7, 0 ), wxGBSpan( 1, 1 ), wxALL|wxEXPAND, 5 );
	
	m_staticText1 = new wxStaticText( this, wxID_ANY, wxT("Geometry Parameters"), wxDefaultPosition, wxDefaultSize, wxALIGN_CENTRE );
	m_staticText1->Wrap( -1 );
	m_staticText1->SetFont( wxFont( 10, 74, 90, 92, false, wxT("Arial") ) );
	
	MainSizer->Add( m_staticText1, wxGBPosition( 0, 0 ), wxGBSpan( 1, 1 ), wxALL, 5 );
	
	
	this->SetSizer( MainSizer );
	this->Layout();
	
	this->Centre( wxBOTH );
}

MainFrame::~MainFrame()
{
}
