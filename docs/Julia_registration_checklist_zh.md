# BoundaryTMatrixMethods.jl Julia Registration Checklist

This note records the practical release path for making the package installable
through Julia's package manager.

## Current Install Modes

Before General registry registration:

```julia
import Pkg
Pkg.add(url = "https://github.com/JuliaRemoteSensing/BoundaryTMatrixMethods.jl")
```

After General registry registration:

```julia
import Pkg
Pkg.add("BoundaryTMatrixMethods")
```

## Required Before General Registration

- Public GitHub repository:
  `https://github.com/JuliaRemoteSensing/BoundaryTMatrixMethods.jl`
- Valid `Project.toml` with stable package name, UUID, version, deps, and compat.
- OSI-compatible license file.
- Passing `Pkg.test()` in a clean environment.
- README with install and quick-start examples.
- CI workflow for Linux, Windows, and macOS.
- All non-standard-library dependencies resolvable from the registries used by
  end users.

The last point is the key blocker to check carefully. If `TransitionMatrices.jl`
is not registered in Julia General, then `BoundaryTMatrixMethods.jl` cannot be
installed by name from General without first registering `TransitionMatrices.jl`
or providing a separate registry that contains both packages.

## Recommended Order

1. Push this package to GitHub and verify URL installation.
2. Confirm whether `TransitionMatrices.jl` is registered in Julia General.
3. If needed, register `TransitionMatrices.jl` first or create a
   `JuliaRemoteSensing` registry containing both packages.
4. Once dependencies are registry-resolvable, install JuliaRegistrator on the
   repository.
5. Trigger registration from a GitHub issue or commit comment:

   ```text
   @JuliaRegistrator register
   ```

6. Wait for the generated PR in `JuliaRegistries/General` to pass AutoMerge.
7. Enable TagBot so a Git tag is created after registration is merged.

## CPC/Software-Paper Positioning

For a CPC-style release, the package should not claim SVM or PMM as new
algorithms. The software contribution is the Julia implementation, the
`TransitionMatrices.jl`-compatible T-matrix interface, the residual diagnostics,
and the conservative PMM/IITM selector that documents accepted and rejected
parameter regimes.
