#include "Cuts.h"
#include "Categories.h"

namespace Selections
{
  const SelectionParams ncpizero_incl = { "ncpizero_incl",
                                          Cuts::ncpizero_incl_cuts,
                                          Cuts::ncpizero_incl_broad_cuts,
                                          Categories::ncpizero_incl_categories,
                                          Categories::true_ncpizero_incl_cut
  };

  const SelectionParams ncpizero_0p0pi = { "ncpizero_0p0pi",
                                           Cuts::ncpizero_0p0pi_cuts,
                                           Cuts::ncpizero_0p0pi_broad_cuts,
                                           Categories::ncpizero_0p0pi_categories,
                                           Categories::true_ncpizero_0p0pi_cut
  };

  const SelectionParams ncpizero_Np0pi = { "ncpizero_Np0pi",
                                           Cuts::ncpizero_Np0pi_cuts,
                                           Cuts::ncpizero_Np0pi_broad_cuts,
                                           Categories::ncpizero_Np0pi_categories,
                                           Categories::true_ncpizero_Np0pi_cut
  };
}
