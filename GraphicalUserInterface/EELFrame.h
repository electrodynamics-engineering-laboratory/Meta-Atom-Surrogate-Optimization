#include <wx/wxprec.h>
#ifndef WX_PRECOMP
 #include <wx/wx.h>
 #include <wx/sizer.h>
#endif

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
