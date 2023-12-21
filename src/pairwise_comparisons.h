#ifndef PAIRWISE_H
#define PAIRWISE_H

arma::vec propose_pairwise_augmentation(
    const arma::vec& ranking,
    const std::vector<std::vector<unsigned int>>& items_above,
    const std::vector<std::vector<unsigned int>> & items_below);

arma::vec propose_swap(
    const arma::vec& ranking,
    const std::vector<std::vector<unsigned int>>& items_above,
    const std::vector<std::vector<unsigned int>>& items_below,
    int& g_diff,
    const int& swap_leap);

#endif
