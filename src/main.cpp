
#include "kernels/graphkernels.h"
#include <Dense>
#include <iostream>

int main() {
  GEDEditCosts edit = GEDEditCosts::MABTS;
  GEDMethods method = GEDMethods::BIPARTITE;
  GEDKernel ged(edit, method);
  // GEDKernel ged();

  return 0;
}