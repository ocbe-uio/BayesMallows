#pragma once

#include "classes.h"

arma::vec propose_pairwise_augmentation(
    const arma::vec& ranking,
    const doubly_nested& items_above,
    const doubly_nested & items_below);

arma::vec propose_swap(
    const arma::vec& ranking,
    const doubly_nested& items_above,
    const doubly_nested& items_below,
    int& g_diff,
    int swap_leap);
