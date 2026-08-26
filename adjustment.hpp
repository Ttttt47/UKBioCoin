#ifndef UKBIOCOIN_ADJUSTMENT_HPP
#define UKBIOCOIN_ADJUSTMENT_HPP

#include <string>

namespace ukc {

struct AdjustmentOptions {
    std::string input_path;
    std::string output_prefix;
    std::string output_format = "table";
    bool lambda_supplied = false;
    double lambda = 1.0;
};

void run_adjustment(const AdjustmentOptions &options);

}  // namespace ukc

#endif
