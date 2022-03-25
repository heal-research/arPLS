#pragma once

#include "linalg.h"
#include "solvers.h"
#include "ap.h"
#include "statistics.h"
#include <math.h>
#include <iostream>
#include <fstream>

alglib::real_1d_array smooth(alglib::real_1d_array& y, double lambda, double ratio);