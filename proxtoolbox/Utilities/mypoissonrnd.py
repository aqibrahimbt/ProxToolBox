from numpy import round, exp, sqrt,zeros_like, ones_like, where, random, ones

def mypoissonrnd(lambd):
    L = exp(-lambd)
    n = zeros_like(lambd)
    k=0
    maxit = 2000
    p = ones_like(lambd)

    ind = where((p>L) & (lambd<200))
    large_ind = where(lambd>=200)

    while (ind[0].size != 0) and k<=maxit:
        print(ind[0].size)
        p[ind] = p[ind] * random.rand(ind[0].size)
        n= n + (p>L)
        ind = where(p>L)
        k += 1

    n[large_ind] = round(lambd[large_ind]+ sqrt(lambd[large_ind])*random.randn(large_ind[0].size))
    neg_ind = where(n<=0)
    n[neg_ind] = 0  

    return n
    