/*
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"

int main(){
  std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();
  matlab::data::ArrayFactory factory;

  
  return 0;
}
*/

#include <wx/wxprec.h>

#ifndef WX_PRECOMP
     #include <wx/wx.h>
     #include <wx/sizer.h>
#endif


//Declare EELImage class for panel
class EELImagePanel: public wxPanel{
  wxBitmap image;
public:
  wxImagePanel(wxFrame* parent, wxString file, wxBitmapType format);

  void paintEvent(wxPaintEvent & evt);
  void paintNow();
  void render(wxDC& dc);

  DECLARE_EVENT_TABLE();
};

BEGIN_EVENT_TABLE(EELImagePanel, wxPanel)
EVT_PAINT(EELImagePanel::paintEvent)
END_EVENT_TABLE();

EELImagePanel::wxImagePanel(wxFrame* parent, wxString file, wxBitmapType format):wxpanel(parent){
  image.LoadFile(file, format);
}

void EELImagePanel::paintEvent(wxPaintEvent & evt){
  wxPaintDC dc(this);
  render(dc);
}

void EELImagePanel::paintNow(){
  wxClientDC dc(this);
  render(dc);
}

void EELImagePanel::render(wxDC& dc){
  dc.DrawBitmap(image, 0, 0, false);
}

class EELApp: public wxApp{
  EELFrame *frame;
  EELImagePanel * drawPane;
public:
  virtual bool OnInit();
};

class EELFrame: public wxFrame{
public:
  EELFrame(const wxString& title, const wxPoint& pos, const wxSize& size);

private:
  void OnHello(wxCommandEvent& event);
  void OnFunctionOne(wxCommandEvent& event);
  void OnFunctionTwo(wxCommandEvent& event);
  void OnExit(wxCommandEvent& event);
  void OnAbout(wxCommandEvent& event);

  wxDECLARE_EVENT_TABLE();
};

wxBEGIN_EVENT_TABLE(EELFrame, wxFrame)
EVT_MENU(1, EELFrame::OnFunctionOne)
EVT_MENU(2, EELFrame::OnFunctionTwo) //Seems to be overwriting previous event menu table value
EVT_MENU(wxID_EXIT, EELFrame::OnExit)
EVT_MENU(wxID_ABOUT, EELFrame::OnAbout)
wxEND_EVENT_TABLE()

wxIMPLEMENT_APP(EELApp);

bool EELApp::OnInit(){

  wxInitAllImageHandlers();

  wxBoxSizer* sizer = new wxBoxSizer(wxHorizontal);
  
  EELFrame *frame = new EELFrame("Electrodynamics Engineering Laboratory Metamodel Suite", wxPoint(0,0), wxSize(2000,1200));

  size->Add(drawPane, 1, wxEXPAND);

  frame->SetSize(sizer);
  frame->Show(true);
  return true;
}

EELFrame::EELFrame(const wxString& title, const wxPoint& pos, const wxSize& size): wxFrame(NULL, wxID_ANY, title, pos, size){
  wxMenu *menuFile = new wxMenu;
  menuFile->AppendSeparator();
 
  menuFile->Append(wxID_EXIT);

  wxMenu *menuHelp = new wxMenu;
  menuHelp->Append(wxID_ABOUT);

  wxMenu *menuFunction = new wxMenu;
  menuFunction->Append(2, "&Function 1\tCtrl-T", "Engage Function 1");
  menuFunction->Append(2, "&Function 2\tCtrl-Y", "Engage Function 2");
   
  wxMenuBar *menuBar = new wxMenuBar;
  menuBar->Append(menuFile, "&File");
  menuBar->Append(menuFunction, "&Function");
  menuBar->Append(menuHelp, "&Help");

  SetMenuBar(menuBar);

  CreateStatusBar();
  SetStatusText("Electrodynamics Engineering Laboratory Metamodel Suite");
}

void EELFrame::OnExit(wxCommandEvent& event){
  Close(true);
}

void EELFrame::OnAbout(wxCommandEvent& event){
  wxMessageBox("This is the EEL Metamodel Suite. See the full documentation for further descriptions of all functions available.", "About Hello World", wxOK | wxICON_INFORMATION);
}

void EELFrame::OnHello(wxCommandEvent& event){
  wxLogMessage("ALL HAIL BOB!");
}

void EELFrame::OnFunctionOne(wxCommandEvent& event){
  wxMessageBox("Function One is not yet implemented","Function One", wxOK | wxICON_INFORMATION);
}

void EELFrame::OnFunctionTwo(wxCommandEvent& event){
  wxMessageBox("Function Two is not yet implemented", "Function Two", wxOK | wxICON_INFORMATION);
}
