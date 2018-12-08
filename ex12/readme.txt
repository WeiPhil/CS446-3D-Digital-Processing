Workload repartition:
    Philippe Weier      34%
    Lucas Ramirez       33%
    Maxime Pisa         33%

----------------------------------

- Pre-processing of our cat model

The first steps we have done on our model were not related to any geometry processing operations, we wanted to 
have a model without any additional spurious artifacts. What we mean by "artifacts" is best shown on the image 
'cat_unpreprocessed.png'. So we removed those artifacts manually using blender giving us the resulting 3D 
object depictured in 'cat_preprocessed.png'.

In the coming geometric processing operations we will always work with our pre-processed model as a basis namely
'cat_highres_preprocessed.obj'

- Curvature Estimation (mean and Gaussian)

