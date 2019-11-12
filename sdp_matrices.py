For cluster 1:
array([[ 1.00000000e+00, -3.98197710e-12,  6.66666670e-01,  3.33333330e-01,  6.66666670e-01,  3.33333330e-01],
       [-3.98197710e-12,  1.30356449e-10,  2.15495973e-11,  2.19060201e-11,  2.15495764e-11,  2.19059989e-11],
       [ 6.66666670e-01,  2.15495973e-11,  5.98720439e-01,  6.79462311e-02,  3.33333340e-01,  3.33333330e-01],
       [ 3.33333330e-01,  2.19060201e-11,  6.79462311e-02,  2.65387099e-01,  3.33333330e-01, -4.52748782e-11],
       [ 6.66666670e-01,  2.15495764e-11,  3.33333340e-01,  3.33333330e-01,  5.98720434e-01,  6.79462356e-02],
       [ 3.33333330e-01,  2.19059989e-11,  3.33333330e-01, -4.52748782e-11,  6.79462356e-02,  2.65387095e-01]])

For clusters 2 & 3:
array([[ 5.98720438e-01,  6.79462320e-02,  6.66666670e-01,  2.15495940e-11,  3.33333340e-01,  3.33333330e-01],
       [ 6.79462320e-02,  2.65387098e-01,  3.33333330e-01,  2.19060168e-11,  3.33333330e-01, -4.52753710e-11],
       [ 6.66666670e-01,  3.33333330e-01,  1.00000000e+00, -3.98197808e-12,  6.66666670e-01,  3.33333330e-01],
       [ 2.15495940e-11,  2.19060168e-11, -3.98197808e-12,  1.30356447e-10,  2.15495841e-11,  2.19060067e-11],
       [ 3.33333340e-01,  3.33333330e-01,  6.66666670e-01,  2.15495841e-11,  5.98720438e-01,  6.79462319e-02],
       [ 3.33333330e-01, -4.52753710e-11,  3.33333330e-01,  2.19060067e-11,  6.79462319e-02,  2.65387098e-01]])
a 1, a-1, a 2, a-2, a 3, a-3

# what do these matrices actually represent? (they're positive semidefinite, but seem weird in form)
# why do they have negative entries?
# rank 2
# every 2x2 submatrices have sum to 1
# imposing that eigenvalues must all be less than 3?
# trace of a 6x6 matrix? it's going to be bounded above


