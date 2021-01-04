data0=sample_mallows(rho0=1:10,alpha=20,n_samples=1000)

rank_distance_vec(rankings=data0,rho=1:10,metric="kendall")
