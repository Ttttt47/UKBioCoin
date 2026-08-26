#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include "statistics.hpp"

namespace {

bool close(double actual, double expected, double tolerance = 2e-12) {
    return std::abs(actual - expected) <=
        tolerance * std::max(1.0, std::abs(expected));
}

int fail(const char *message, double actual = 0.0) {
    std::cerr << message;
    if (actual != 0.0) std::cerr << ": " << actual;
    std::cerr << '\n';
    return 1;
}

}  // namespace

int main() {
    if (ukc::residual_degrees_of_freedom(100, 0) != 98 ||
        ukc::residual_degrees_of_freedom(71, 1) != 68 ||
        ukc::residual_degrees_of_freedom(4, 1) != 1 ||
        ukc::residual_degrees_of_freedom(3, 1) != 0) {
        return fail("residual df does not equal n-c-2");
    }

    const double t_probability =
        ukc::two_sided_t_neg_log10(2.9192017967990469, 98);
    if (!close(t_probability, 2.3611488418022879, 5e-13)) {
        return fail("two-sided t log probability differs from reference",
                    t_probability);
    }
    if (!std::isnan(ukc::two_sided_t_neg_log10(1.0, 0)) ||
        !std::isnan(ukc::two_sided_t_neg_log10(
            std::numeric_limits<double>::quiet_NaN(), 10))) {
        return fail("invalid t inputs did not return NaN");
    }

    const std::vector<double> z{10.0, 20.0, 40.0, 50.0, 100.0};
    const std::vector<double> expected{
        22.8170234098221,
        88.2590653474116,
        349.135976463682,
        544.665305866333,
        2173.57051287337};
    for (std::size_t i = 0; i < z.size(); ++i) {
        const double actual = ukc::chi_square_one_df_neg_log10_from_z(z[i]);
        if (!close(actual, expected[i], 5e-13)) {
            return fail("chi-square log tail differs from extreme reference",
                        actual);
        }
    }
    if (!close(ukc::chi_square_one_df_neg_log10_from_z(0.0), 0.0) ||
        !std::isnan(ukc::chi_square_one_df_neg_log10_from_z(-1.0))) {
        return fail("chi-square z-domain validation failed");
    }
    return 0;
}
