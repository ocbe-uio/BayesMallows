#ifndef SUBSET_H
#define SUBSET_H

# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>
# include <Rcpp.h>

using namespace std;


int i4_sign ( int i );
void i4vec_decrement ( int n, int v[] );
void perm_ascend ( int n, int a[], int &length, int sub[] );
bool perm0_check ( int n, int p[] );
int perm0_distance ( int n, int a[], int b[] );
int *perm0_inverse ( int n, int p[] );
void perm0_mul ( int n, int p1[], int p2[], int p3[] );


#endif
