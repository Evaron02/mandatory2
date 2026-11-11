# matmek4270-mandatory2

[![MAT-MEK4270 mandatory 2](https://github.com/Evaron02/mandatory2/actions/workflows/main.yml/badge.svg)](https://github.com/Evaron02/mandatory2/actions/workflows/main.yml)

# Some further instructions for this repo

This repository contains a small framework in `galerkin.py` for assembling mass/stiffness-like operators and projecting/solving with different polynomial and trigonometric bases and their boundary-condition lifts. Classes are used for code reuse and clarity, kept as simple as possible for educational purposes.

## How to pass this assignment

There is no manual grading. The assignment is passed once your fork of this repository shows a **green badge**, i.e. `pytest galerkin.py` runs `test_project`, `test_convection_diffusion` and `test_helmholtz` without errors on GitHub Actions.

Check locally before pushing:
```bash
uv sync
uv run ty check galerkin.py   # optional, but should pass — see "Type hints" below
uv run pytest -v galerkin.py
# or simply:
uv run python galerkin.py
```

## What you need to implement

`galerkin.py` marks unfinished code two different ways, and they mean different things:

- **`NotImplementedError`** marks a concrete method for you to fill in — mirror the pattern of the sibling classes/methods that are already implemented.
- **`@abstractmethod`** (from `abc`) marks a method that's never meant to be implemented at the class where it's declared — only overridden by a concrete subclass further down the hierarchy, which every subclass eventually does (directly, or via one of the `NotImplementedError` blanks below). If you see `@abstractmethod` above a method with a bare `...` body, **leave it alone**.

Implement:
- [X] `Polynomial.reference_domain` — replaced Not Implimented Error by the domain (-1,1), seems to work just fine :)
- [X] `Polynomial.derivative_basis_function` — Implemented the k'th derivative of `basis_function(j)`; did a trial run using 
  ```bash
  uv run python -c "from galerkin import Legendre; print(Legendre(4).derivative_basis_function(2, 2)(0))"
  ```
    which printed 
    ```text
    3.0 
    ```
    so it seems to be working :D
- [X] `Polynomial.mass_matrix` — Implemented the (diagonal) mass matrix built from `self.L2_norm_sq`, should work for `Legendre` and `Chebyshev`. did a trial run by using 
  ```bash
  uv run python -c "from galerkin import Legendre; V = Legendre(3); V.L2_norm_sq = lambda n: [2.0, 4.0, 6.0]; print(V.mass_matrix().toarray())"
  ```
    which printed
  ```text
  [[2. 0. 0.]
   [0. 4. 0.]
   [0. 0. 6.]]
  ```
  Which is correct, so this method should be working jsut fine :)
- [X] `Legendre.L2_norm_sq` and `Chebyshev.L2_norm_sq` — Implemented both squared norms, according to they're properties. to check if it worked i did:
  ```bash
  uv run python -c "from galerkin import Legendre, Chebyshev; print(Legendre(5).L2_norm_sq(3)); print(Chebyshev(5).L2_norm_sq(3))"
  ```
    which printed
  ```text
  [2.         0.66666667 0.4       ]
  [3.14159265 1.57079633 1.57079633]
  ```
  seems correct, so these two methods should be correct.
  I added Bounds for the `Chebyshev.L2_norm_sq` method, because it needed to handle case `N = 0` and `N > 0`, So i added one for the other method too, but these shouldne be nescessary 
- [ ] `Chebyshev.eval` — mirror `Legendre.eval`, using `np.polynomial.chebyshev.chebval` instead of `legval`z.
- [ ] `Cosines` (`__init__`, `basis_function`, `derivative_basis_function`, `L2_norm_sq`) — mirror `Sines`, but with `cos(jπx)` and its `j=0` mode handled separately (its derivative vanishes identically for any order)z.
- [ ] `DirichletLegendre.basis_function` — mirror `DirichletChebyshev.basis_function`, using Legendre polynomials instead of Chebyshevz.
- [ ] `NeumannLegendre` and `NeumannChebyshev` (`__init__`, `basis_function`) — composite spaces enforcing homogeneous derivative (Neumann) boundary conditions via a stencil matrix `S`, following the same pattern as `DirichletLegendre`/`DirichletChebyshev`, but solving for a stencil coefficient that makes `ψ_j'(±1) = 0`. See [lecture 11](https://matmek-4270.github.io/matmek4270-book/lecture11.html#neumann-boundary-conditions) for the boundary-derivative values to use.

## Architecture, in brief

The module docstring and class docstrings in `galerkin.py` describe the architecture — mapping helpers, `FunctionSpace` and its orthogonal/trigonometric/composite subclasses, the `Dirichlet`/`Neumann` boundary lifts, and the `TrialFunction`/`TestFunction`/`inner` weak-form domain-specific language (DSL). Start there (`python -c "help(__import__('galerkin'))"`, or just read top-to-bottom) before diving into the blanks above.

## Gotchas

- **Domain mapping.** Always go through `map_reference_domain`/`map_true_domain` before evaluating a basis function or its derivative — don't evaluate directly in the physical domain.
- **Chebyshev weight.** `Chebyshev.inner_product` substitutes `x = cos(θ)` and integrates over `[0, π]` to handle the `1/√(1−x²)` weight; don't use plain `quad` on `[-1, 1]` there. `assemble_generic_matrix` instead passes `weight="alg", wvar=(-0.5, -0.5)` to `quad` for the same reason.
- **Coefficient lengths.** Simple spaces use `N+1` coefficients; composite spaces use a stencil `S` of shape `(N+1, N+3)` and a diagonal `M` of length `N+3` in the underlying orthogonal basis.
- **Trigonometric derivatives.** `Sines.derivative_basis_function` shows the scale/sign pattern for `((j+1)π)^k`; `Cosines` needs the analogous pattern for `(jπ)^k`, with the `j=0` mode's derivative identically zero for any order.
- **Boundary lifts in the RHS.** When assembling `b = inner(f - operator(B.x), v)`, remember to subtract the same differential operator applied to the lift `B.x`, as shown in `test_helmholtz`/`test_convection_diffusion`.
- **Type hints & `@overload`.** Method signatures are fully type-annotated, so they describe exactly what each blank should accept and return. `basis_function(j, sympy=False)` returns a different type depending on `sympy` — a callable when `sympy=False`, a sympy expression when `sympy=True` — which a plain type hint can't express precisely. It's declared with two `@overload` stubs (one per literal value of `sympy`), followed by the actual (untyped) implementation; only the two stubs are seen by the type checker, so it can pick the right return type at each call site instead of an imprecise union of both. Your implementation should keep `uv run ty check galerkin.py` passing — it catches a wrong return type or shape immediately, well before the (much slower) numerical tests would.
