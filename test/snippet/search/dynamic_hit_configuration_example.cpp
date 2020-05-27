#include <seqan3/search/configuration/all.hpp>

int main()
{
    // Default constructed: Has no hit strategy selected.
    seqan3::search_cfg::hit dynamic_hit{};

    // Select hit_all
    dynamic_hit = seqan3::search_cfg::hit_all;

    // If condition is true chose strata strategy, otherwise find the single best hit.
    if (true)
        dynamic_hit = seqan3::search_cfg::hit_strata{4};
    else
        dynamic_hit = seqan3::search_cfg::hit_single_best;

    // Combine it with other configurations.
    seqan3::configuration const cfg = dynamic_hit | seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}};

    // Directly initialised.
    seqan3::search_cfg::hit dynamic_hit2{seqan3::search_cfg::hit_all_best};

    // You cannot combine the dynamic hit configuration with the static ones.
    // auto fail = seqan3::search_cfg::hit_single_best | seqan3::search_cfg::hit; // doesn't compile

    return 0;
}
