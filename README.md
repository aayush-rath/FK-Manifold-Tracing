# Manifold Tracing using FK-Triangulation

This repository contains an independent implementation of the Freudenthal–Kuhn triangulation–based manifold tracing algorithm described in the PhD thesis “Meshing Submanifolds using Coxeter Triangulation” by Siargey Kachanovich.

The implementation focuses specifically on tracing co-dimension 1 manifolds, i.e. manifolds of dimension 
$d−1$ embedded in an ambient space of dimension $d$, using a permutahedral representation of the Freudenthal–Kuhn triangulation.

The algorithmic formulation follows the presentation in Kachanovich’s thesis and is conceptually similar to the approach used in the GUDHI library (MIT licensed), but the code in this repository is written from scratch and does not depend on GUDHI.

For a mature and production-ready implementation of the same algorithmic framework, please refer to the
[GUDHI documentation](https://gudhi.inria.fr/doc/latest/)
