% Number of threads = 1

% number of repeats:% 3
% enter first, last, inc:% 48 1968 48 
data = [
%  n          reference      |         current implementation 
%        time       GFLOPS   |    time       GFLOPS     diff 
  1968 3.0676e-01 4.9694e+01    3.0410e-01 5.0129e+01 3.0127e-12
  1920 2.8534e-01 4.9609e+01    2.8679e-01 4.9359e+01 2.9559e-12
  1872 2.6581e-01 4.9360e+01    2.6274e-01 4.9936e+01 2.7285e-12
  1824 2.4723e-01 4.9092e+01    2.4393e-01 4.9756e+01 2.7285e-12
  1776 2.2585e-01 4.9606e+01    2.2345e-01 5.0140e+01 2.6148e-12
  1728 2.0830e-01 4.9543e+01    2.0708e-01 4.9833e+01 2.6716e-12
  1680 1.9393e-01 4.8900e+01    1.9121e-01 4.9597e+01 2.4443e-12
  1632 1.7750e-01 4.8977e+01    1.7502e-01 4.9672e+01 2.5011e-12
  1584 1.6322e-01 4.8699e+01    1.6033e-01 4.9578e+01 2.1600e-12
  1536 1.4572e-01 4.9736e+01    1.4527e-01 4.9893e+01 2.2169e-12
  1488 1.3393e-01 4.9199e+01    1.3317e-01 4.9481e+01 2.1600e-12
  1440 1.2209e-01 4.8914e+01    1.2129e-01 4.9235e+01 1.9895e-12
  1392 1.1043e-01 4.8848e+01    1.0973e-01 4.9161e+01 1.8758e-12
  1344 1.0061e-01 4.8262e+01    9.8954e-02 4.9068e+01 1.8190e-12
  1296 9.1097e-02 4.7790e+01    8.8283e-02 4.9314e+01 1.7621e-12
  1248 7.9561e-02 4.8862e+01    7.9340e-02 4.8998e+01 1.5916e-12
  1200 7.1082e-02 4.8620e+01    7.0153e-02 4.9264e+01 1.6485e-12
  1152 6.2377e-02 4.9019e+01    6.1933e-02 4.9371e+01 1.5348e-12
  1104 5.5736e-02 4.8284e+01    5.5114e-02 4.8829e+01 1.2506e-12
  1056 4.8743e-02 4.8318e+01    4.7995e-02 4.9071e+01 1.2506e-12
  1008 4.1468e-02 4.9397e+01    4.1297e-02 4.9602e+01 1.1937e-12
   960 3.6123e-02 4.8984e+01    3.5860e-02 4.9344e+01 1.1084e-12
   912 3.1137e-02 4.8723e+01    3.0808e-02 4.9243e+01 9.3792e-13
   864 2.6777e-02 4.8173e+01    2.6358e-02 4.8939e+01 8.8107e-13
   816 2.2620e-02 4.8040e+01    2.2154e-02 4.9051e+01 7.6739e-13
   768 1.8107e-02 5.0035e+01    1.8381e-02 4.9288e+01 7.9581e-13
   720 1.5316e-02 4.8740e+01    1.5023e-02 4.9689e+01 6.5370e-13
   672 1.2320e-02 4.9262e+01    1.2210e-02 4.9706e+01 6.5370e-13
   624 9.8080e-03 4.9546e+01    9.7912e-03 4.9630e+01 5.6843e-13
   576 7.6730e-03 4.9812e+01    7.6657e-03 4.9859e+01 5.1159e-13
   528 5.9398e-03 4.9563e+01    5.8615e-03 5.0225e+01 3.6948e-13
   480 4.4332e-03 4.9893e+01    4.4304e-03 4.9924e+01 2.9843e-13
   432 3.2287e-03 4.9941e+01    3.2243e-03 5.0008e+01 2.7001e-13
   384 2.3025e-03 4.9185e+01    2.2725e-03 4.9834e+01 2.1316e-13
   336 1.5595e-03 4.8647e+01    1.5345e-03 4.9440e+01 1.7053e-13
   288 9.9393e-04 4.8067e+01    9.7129e-04 4.9188e+01 1.1369e-13
   240 5.7269e-04 4.8277e+01    5.6616e-04 4.8834e+01 4.2633e-14
   192 3.0183e-04 4.6900e+01    2.9623e-04 4.7787e+01 2.8422e-14
   144 1.3334e-04 4.4788e+01    1.2821e-04 4.6579e+01 2.8422e-14
    96 4.4062e-05 4.0159e+01    3.9945e-05 4.4298e+01 1.0658e-14
    48 9.0170e-06 2.4530e+01    5.9010e-06 3.7482e+01 5.3291e-15
];

% Maximum difference between reference and your implementation: 3.012701e-12.
