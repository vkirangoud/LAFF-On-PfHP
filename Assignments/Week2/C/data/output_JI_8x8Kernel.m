% number of repeats:% 1000
% enter first, last, inc:% 48 1488 48 
data = [
%  n          reference      |         current implementation 
%        time       GFLOPS   |    time       GFLOPS     diff 
  1488 1.3357e-01 4.9332e+01    4.7862e-01 1.3767e+01 2.1600e-12
  1440 1.2228e-01 4.8838e+01    2.0015e-01 2.9838e+01 2.1032e-12
  1392 1.1051e-01 4.8815e+01    1.8399e-01 2.9320e+01 1.9327e-12
  1344 1.0016e-01 4.8475e+01    1.6303e-01 2.9783e+01 1.7621e-12
  1296 9.0098e-02 4.8321e+01    1.4393e-01 3.0249e+01 1.8190e-12
  1248 7.8863e-02 4.9295e+01    1.3024e-01 2.9849e+01 1.6485e-12
  1200 7.0505e-02 4.9018e+01    1.1494e-01 3.0067e+01 1.5348e-12
  1152 6.2840e-02 4.8658e+01    1.0152e-01 3.0118e+01 1.3642e-12
  1104 5.5169e-02 4.8780e+01    8.9539e-02 3.0055e+01 1.2506e-12
  1056 4.8476e-02 4.8584e+01    7.7583e-02 3.0357e+01 1.1937e-12
  1008 4.1639e-02 4.9194e+01    6.7852e-02 3.0189e+01 1.0232e-12
   960 3.6357e-02 4.8669e+01    5.8611e-02 3.0190e+01 1.0232e-12
   912 3.1095e-02 4.8788e+01    4.9896e-02 3.0405e+01 1.0800e-12
   864 2.6300e-02 4.9047e+01    4.2384e-02 3.0435e+01 8.5265e-13
   816 2.1940e-02 4.9530e+01    3.5703e-02 3.0437e+01 8.2423e-13
   768 1.8033e-02 5.0238e+01    3.0192e-02 3.0007e+01 7.1054e-13
   720 1.4917e-02 5.0043e+01    2.4494e-02 3.0476e+01 6.8212e-13
   672 1.2068e-02 5.0291e+01    1.9918e-02 3.0471e+01 5.6843e-13
   624 9.6522e-03 5.0345e+01    1.5942e-02 3.0483e+01 5.4001e-13
   576 7.6436e-03 5.0003e+01    1.2577e-02 3.0390e+01 4.8317e-13
   528 5.9305e-03 4.9641e+01    9.6771e-03 3.0422e+01 4.2633e-13
   480 4.4010e-03 5.0258e+01    7.2740e-03 3.0407e+01 3.1264e-13
   432 3.2241e-03 5.0012e+01    5.3033e-03 3.0404e+01 2.7001e-13
   384 3.3503e-03 3.3802e+01    3.7246e-03 3.0405e+01 5.6843e-14
   336 2.2518e-03 3.3692e+01    2.4889e-03 3.0482e+01 5.6843e-14
   288 1.4234e-03 3.3566e+01    1.5677e-03 3.0476e+01 7.1054e-14
   240 5.4896e-04 5.0364e+01    9.0926e-04 3.0407e+01 4.2633e-14
   192 2.8006e-04 5.0546e+01    4.6249e-04 3.0608e+01 2.8422e-14
   144 1.1636e-04 5.1324e+01    1.9548e-04 3.0551e+01 2.8422e-14
    96 3.5366e-05 5.0033e+01    5.8199e-05 3.0404e+01 1.0658e-14
    48 4.8780e-06 4.5343e+01    7.2940e-06 3.0324e+01 5.3291e-15
];

% Maximum difference between reference and your implementation: 2.160050e-12.
