#include <wx/wx.h>
#include <wx/sizer.h>

class EELImagePanel: public wxPanel{
  wxImage image;
  wxBitmap resized;
  int w, h;
public:
  EELImagePanel(wxFrame* parent, wxString file, wxBitmapType format);

  void paintEvent(wxPaintEvent& event);
  void paintNow();
  void OnSize(wxSizeEvent& event);
  void render(wxDC& dc);

  DECLARE_EVENT_TABLE()
  
};

BEGIN_EVENT_TABLE(EELImagePanel, wxPanel)
EVT_PAINT(EELImagePanel::paintEvent) //Catch paint events
EVT_SIZE(EELImagePanel::OnSize) //Size event
END_EVENT_TABLE()
