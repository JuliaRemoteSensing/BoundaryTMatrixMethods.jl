# Contributing

Contributions are welcome through issues and pull requests.

For numerical changes, please include a small validation case or explain which
existing SVM, PMM, EBCM, or IITM comparison is affected. Changes that alter
T-matrix construction should keep the public return type compatible with
`TransitionMatrices.jl`.

Before opening a pull request, run:

```julia
import Pkg
Pkg.test()
```
