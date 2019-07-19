/*
#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"

int main(){
  std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();
  matlab::data::ArrayFactory factory;

  
  return 0;
}
*/

#include "ElectrodynamicsEngineeringLaboratoryGUI.h"
#include "EELFrame.h"


bool EELApp::OnInit(){

  wxInitAllImageHandlers();
  
  EELFrame *frame = new EELFrame("Electrodynamics Engineering Laboratory Metamodel Suite", wxPoint(0,0), wxSize(2000,1200));

  frame->Show(true);
  return true;
}
