# ArnoldiMethodConvenienceWrapper.jl

Convenience wrapper for ArnoldiMethod.jl. 

- The function `eigs` that is exported takes the same arguments as the eponymous function from Arpack.jl.
- Only the calculation of the lowest-magnitude eigenvalue is supported.
