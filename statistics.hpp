#ifndef UKBIOCOIN_STATISTICS_HPP
#define UKBIOCOIN_STATISTICS_HPP

#include <cstddef>

namespace ukc {

int residual_degrees_of_freedom(int sample_size,
                                std::size_t covariate_count);

double two_sided_t_neg_log10(double t_statistic, int degrees_of_freedom);

double chi_square_one_df_neg_log10(double chi_square);

double chi_square_one_df_neg_log10_from_z(double absolute_z);

}  // namespace ukc

#endif
