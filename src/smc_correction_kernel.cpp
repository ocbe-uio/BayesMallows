
  //' @description Function to determine if the augmented ranking is compatible with the new observed partial ranking.
  //'              If it is not, the we create a new augmentation using the random sampling approachand calculate the augmentation probability.
  //'
  //' INPUT:
  //' @param current_ranking A ranking sequence vector of the current augmented ranking (no missing values)
  //' @param observed_ranking  A ranking sequence vector of the observed partial ranking (no missing values) The original incomplete partial ranking is in the rankings data set.
  //' @param n_items Integer is the number of items in a ranking

  //' OUTPUT: List containing the proposed 'corrected' augmented ranking that is compatible with the new observed ranking for a user


Rcpp::List (int n_items,
    arma::vec observed_ranking, //R_obs
    arma::vec current_ranking) {  // R_curr

          //------------------------------------------------------------------------------------------
          //TO FIX
          check if new information means 'mistakes' made with augmented rankings
          check = (observed_ranking == current_ranking);
          bool condition_1 = any(check == FALSE, na.rm = TRUE) == FALSE;
          if (condition_1) {
          //------------------------------------------------------------------------------------------

             double correction_prob = 1.0

                 //return(output)
                 return Rcpp::List::create(
                     Rcpp::Named("ranking") = current_ranking,
                     Rcpp::Named("correction_prob") = correction_prob
                 );


         }
         else {

             // resample from smaller pool of possible augmented rankings
             // select uniform the proposed ranking compatible with the known observed rankings
             //ranks = c(1:n_items)

             //  find elements missing from original observed ranking
             arma::vec remaining_set = current_ranking.elem(find_nonfinite(observed_ranking));
             //remaining_set = unique(ranks[!ranks % in % R_obs])

             // create new agumented ranking by sampling remaining ranks from set uniformly
             arma::vec proposed_ranking = observed_ranking;

             if (remaining_set.n_elem == 1) {
                 //R_prop[is.na(R_prop)] = remaining_set
                 arma::uvec unranked_items = find_nonfinite(proposed_ranking);
                 proposed_ranking.elem(unranked_items) = remaining_set;

             }
             else {
                 //R_prop[is.na(R_prop)] <- sample(remaining_set, size = sum(is.na(R_prop)), replace = F)
                 // generate random order for remaining_set
                 remaining_set = std::move(arma::shuffle(remaining_set));
                 proposed_ranking.elem(unranked_items) = remaining_set;

             };

             double correction_prob = 1.0 / factorial(remaining_set.n_elem);

             //return(output)
             return Rcpp::List::create(
                 Rcpp::Named("ranking") = proposed_ranking,
                 Rcpp::Named("correction_prob") = correction_prob
             );


         }

 }
