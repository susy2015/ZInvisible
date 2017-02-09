#include "ScaleFactorsttBar.h"

extern "C"
{
    static double python_sf_norm0b()
    {
        return ScaleFactors::sf_norm0b();
    }
    static double python_sf_norm0b_err()
    {
        return ScaleFactors::sfunc_norm0b();
    }
}
