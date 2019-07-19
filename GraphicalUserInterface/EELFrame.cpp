#include "EELFrame.h"

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
