Workload repartition:
    Philippe Weier      34%
    Lucas Ramirez       33%
    Maxime Pisa         33%

----------------------------------

- Comparing Lx = b and L2x = b

We see that for the non-squared laplacian the model behaves more like 
a minimal surface in the sense that it tries to minimize the deformation energy.
The model streches very locally.

With the Laplacian squared, the whole model bends because it tries 
to resist more like thin plates by not streching but more by bending.

- Comparing Different Laplacian Weights

The uniform weights do not really respect the original curvature once the bending is done.
Using these weights seems to result in a bigger loss of information compared to cotangent weights. 

- Comparing with Physical Deformation

The difference comes when we deform to the extreme: real materials would snap or break whereas 
here any deformation is computed without the object "breaking". Another difference is self-intresection, 
such things aren't possible in the physical world but they can happen with the methods we use.

----------------------------------

See screenshots in the imgs/ folder