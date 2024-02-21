#pragma once

#include "classes.h"

RankProposal propose_swap(
    const arma::vec& current_rank,
    const doubly_nested& items_above,
    const doubly_nested& items_below,
    int swap_leap);
