class EELApp: public wxApp{
  EELFrame *frame;
  
public:
  virtual bool OnInit();
};

wxBEGIN_EVENT_TABLE(EELFrame, wxFrame)
EVT_MENU(1, EELFrame::OnFunctionOne)
EVT_MENU(2, EELFrame::OnFunctionTwo) //Seems to be overwriting previous event menu table value
EVT_MENU(wxID_EXIT, EELFrame::OnExit)
EVT_MENU(wxID_ABOUT, EELFrame::OnAbout)
wxEND_EVENT_TABLE()

wxIMPLEMENT_APP(EELApp);
