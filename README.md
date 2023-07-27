## GenHurst

### The Hurst exponent

Calculates the generalized Hurst exponent of a time series. The Hurst exponent gives a value indicating the long-term memory of a time-series, similar to the decay of a autocorrelation function: https://en.wikipedia.org/wiki/Hurst_exponent

See also in the code for some further references.

### The implementation

This implementation is a more or less literal translation from Matlab to Python (3.6) of the following Matlab code written by Tomaso Aste in 2013: https://www.mathworks.com/matlabcentral/fileexchange/30076-generalized-hurst-exponent

### How to use the code

A function `mH = genhurst(S,q)` is defined, with `S` the time series to be analyzed as a numpy array and `q` the Hurst exponent to be used, yielding a numerical (mean) value `mH`.

