Workload repartition:
    Philippe Weier      34%
    Lucas Ramirez       33%
    Maxime Pisa         33%

----------------------------------

- 2 Minimal Surfaces

Cylinder 3 differs from the other 2 when we apply our minimal surface
solver to it. After a few iterations the mesh converges to two disks (where
the boundaries originally were) linked by an infinitly thin line.

We observe this result because the cylinder's height is much larger than
its bases diameter, so the "shortest" path to link the 2 bases (where all the
mesh collapses) is an infinitly thin line. This result is consistent with the
goal of minimal surface optimization.

If we replace the cotan Laplacian with the uniform Laplacian the 3 cylinders'
interior vertices "collapse" in a single point. This is because the uniform
Laplacian doesn't take into account geometry but only connectivity. All of the
cylinder's interior (non-boundary) vertices have the same connectivity so they 
collapse to a single (infinitly small) point while the boundary vertices stay
in their original place.

----------------------------------

See screenshots in the imgs/ folder