#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"

int main(){
  std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();
  matlab::data::ArrayFactory factory;

  
  return 0;
}
