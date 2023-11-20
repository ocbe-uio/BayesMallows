# compute posterior intervals works

    Code
      compute_posterior_intervals(m, burnin = 7, level = 0.05, parameter = "alpha")
    Output
        parameter  mean median conf_level          hpdi central_interval
      1     alpha 0.764   0.79        5 % [0.700,0.700]    [0.786,0.791]

---

    Code
      compute_posterior_intervals(m, burnin = 7, level = 0.1, parameter = "cluster_probs")
    Output
            parameter mean median conf_level          hpdi central_interval
      1 cluster_probs    1      1       10 % [1.000,1.000]          [1.000]

---

    Code
      compute_posterior_intervals(m, burnin = 7, level = 0.01, parameter = "rho")
    Output
         item parameter mean median conf_level hpdi central_interval
      1    P1       rho   17     17        1 % [17]             [17]
      2    P2       rho    7      7        1 %  [7]            [7,7]
      3    P3       rho   16     16        1 % [16]             [16]
      4    P4       rho   15     15        1 % [15]             [15]
      5    P5       rho    2      2        1 %  [2]            [2,2]
      6    P6       rho   19     19        1 % [19]             [19]
      7    P7       rho    1      1        1 %  [1]              [1]
      8    P8       rho   12     12        1 % [12]             [12]
      9    P9       rho   10     10        1 % [10]             [10]
      10  P10       rho    6      6        1 %  [6]            [6,6]
      11  P11       rho    9      9        1 %  [9]              [9]
      12  P12       rho   14     14        1 % [14]             [14]
      13  P13       rho   13     13        1 % [13]             [13]
      14  P14       rho    3      3        1 %  [3]            [3,3]
      15  P15       rho    4      4        1 %  [4]              [4]
      16  P16       rho    8      8        1 %  [8]              [8]
      17  P17       rho   18     18        1 % [18]             [18]
      18  P18       rho    5      5        1 %  [5]              [5]
      19  P19       rho   20     20        1 % [20]             [20]
      20  P20       rho   11     11        1 % [11]             [11]

---

    Code
      compute_posterior_intervals(m, burnin = 8)
    Output
          cluster parameter  mean median conf_level          hpdi central_interval
      1 Cluster 1     alpha 0.666  0.666       95 % [0.621,0.711]    [0.624,0.709]
      2 Cluster 2     alpha 1.122  1.122       95 % [1.065,1.180]    [1.068,1.177]

