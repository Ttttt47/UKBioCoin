#include <cmath>
#include <iostream>

#include "tools/cpp_cdfs-master/cdf_chisqt/cdf_base.h"

int main() {
    const double value = quantile_t(0.75, 1.0);
    if (!std::isfinite(value) || std::abs(value - 1.0) > 1e-10) {
        std::cerr << "quantile_t(0.75, 1) returned " << value << '\n';
        return 1;
    }
    return 0;
}
