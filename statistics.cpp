#include "statistics.hpp"

#include <cmath>
#include <limits>

#include "tools/cpp_cdfs-master/cdf_chisqt/cdf_base.h"

namespace ukc {

int residual_degrees_of_freedom(int sample_size,
                                std::size_t covariate_count) {
    if (covariate_count >
        static_cast<std::size_t>((std::numeric_limits<int>::max)() - 2)) {
        return -1;
    }
    return sample_size - static_cast<int>(covariate_count) - 2;
}

double two_sided_t_neg_log10(double t_statistic, int degrees_of_freedom) {
    if (!std::isfinite(t_statistic) || degrees_of_freedom <= 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double log_p =
        cdf_t_log(-std::abs(t_statistic), degrees_of_freedom) + std::log(2.0);
    return -log_p / std::log(10.0);
}

double chi_square_one_df_neg_log10(double chi_square) {
    if (!std::isfinite(chi_square) || chi_square < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double log_p = cdf_chisq_log(chi_square, 1.0, false);
    return -log_p / std::log(10.0);
}

double chi_square_one_df_neg_log10_from_z(double absolute_z) {
    if (!std::isfinite(absolute_z) || absolute_z < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (absolute_z < 8.0) {
        return chi_square_one_df_neg_log10(absolute_z * absolute_z);
    }

    // With one degree of freedom, P(ChiSq >= z^2) = 2 * Phi(-z).
    // Evaluate the normal-tail expansion directly so the probability never
    // underflows, and scale z before squaring to delay numeric overflow.
    const double inverse_square = 1.0 / (absolute_z * absolute_z);
    double term = 1.0;
    double series = 1.0;
    double previous = std::abs(term);
    for (int order = 1; order < 256; ++order) {
        term *= -(2.0 * order - 1.0) * inverse_square;
        const double magnitude = std::abs(term);
        if (magnitude >= previous) break;
        series += term;
        if (magnitude <= std::numeric_limits<double>::epsilon() *
                             std::abs(series)) {
            break;
        }
        previous = magnitude;
    }
    const double log_10 = std::log(10.0);
    const double scaled = absolute_z / std::sqrt(2.0 * log_10);
    return scaled * scaled +
           (std::log(absolute_z) + 0.5 * std::log(2.0 * std::acos(-1.0)) -
            std::log(series) - std::log(2.0)) /
               log_10;
}

}  // namespace ukc
