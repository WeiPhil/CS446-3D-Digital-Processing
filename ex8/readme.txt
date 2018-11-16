Workload Repartition:
    Maxime Pisa         33%
    Philippe Weier      33%
    Lucas Ramirez       34%


TANGENTIAL SMOOTHING

The algorithm doesn't lead to a stationary solution after a finite
number of steps, i.e. the mesh keeps changing iteration after iteration.

Looking closer at the mesh we see that one of the main reason for change
is the equalize_valences() function. On each remesh some edges are getting
flipped because they keep increasing their endpoints valence in their local
neighborhood.

The tangential_relaxation() function is also a cause of this non-stationarity.
Because the mesh's discrete mean curvature is never 0 everywhere (and not coli-
near to the vertices normals) the vertices position keeps being "smoothed" in
the mesh's tangent plane, creating new opportunities for edge split/collapse.