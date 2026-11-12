# matmek4270-mandatory2 Eva Rønningsbakk

[![MAT-MEK4270 mandatory 2](https://github.com/Evaron02/mandatory2/actions/workflows/main.yml/badge.svg)](https://github.com/Evaron02/mandatory2/actions/workflows/main.yml)

# Everything i did

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
- [X] `Chebyshev.eval` — mirrored `Legendre.eval`, replacing `legval`z with `np.polynomial.chebyshev.chebval` did a little test to see if it works:
  
  ´´´bash
  uv run python -c "import numpy as np; from galerkin import Chebyshev; print(Chebyshev(3, domain=(0, 10)).eval(np.array([1., 2., 3.]), np.array([0., 5., 10.])))"
  ´´´
  which printed
  ``` text
  [ 2. -2.  6.]
  ```
  which is correct, so it seems to be working fine.
- [x] `Cosines` (`__init__`, `basis_function`, `derivative_basis_function`, `L2_norm_sq`) — Implemented the `cos(jπx)` basis, the constant mode has squared norm `1` and zero derivatives of positive order. The remaining modes have squared norm `1/2`.

  Quick check for the squared norms and the constant mode’s first derivative.

  ```bash
  uv run python -c "from galerkin import Cosines; V = Cosines(3); print(V.L2_norm_sq(3)); print(V.derivative_basis_function(0, 1)(0.5) == 0)"
  ```

  this is what printed:

  ```text
  [1.  0.5 0.5]
  True
  ```
  all good here

- [x] `DirichletLegendre.basis_function` — Implemented `P_j - P_(j+2)`, which vanishes at both endpoints of the reference interval.

  checking using this

  ```bash
  uv run python -c "import numpy as np; from galerkin import DirichletLegendre; V = DirichletLegendre(3); print(np.allclose(V.basis_function(1)(np.array([-1., 1.])), 0))"
  ```

  print:

  ```text
  True
  ```
  all good

- [x] `NeumannLegendre` and `NeumannChebyshev` (`__init__`, `basis_function`) — Implemented the stencil matrix `S` and the basis functions so that `ψ_j'(±1) = 0`. 

  quick check that it sseems fine

  ```bash
  uv run python -c "import numpy as np; from galerkin import NeumannLegendre, NeumannChebyshev; points = np.array([-1., 1.]); print(np.allclose(NeumannLegendre(3).derivative_basis_function(1)(points), 0)); print(np.allclose(NeumannChebyshev(3).derivative_basis_function(1)(points), 0))"
  ```

  here is what printed:

  ```text
  True
  True
  ```

### Full test run

```bash
uv run pytest -q galerkin.py  
```
this is what pronted
```text
3 passed, 1 warning in 21.59s
```

The collection warning is for `TestFunction`, a Galerkin helper class whose name starts with `Test`. All three assignment tests passed.
