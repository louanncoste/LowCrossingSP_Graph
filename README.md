A Rust implementation of two algorithms, LowMax and LowSum

# LowMax 

get in entry a graph as a list of edges. The vertices must be numbered from 0 to the size of the graph minus one.
ouputs a path with low maximum crossing number 

parameters: the graph, then the name of the output files then a parameter alpha between 0 and 1 used for the sampling (close to 0 is fast but the path is not very good, close to 1 is long but has a better path)
example:

``` cargo run data/facebook_combined.txt output_alpha_0p5 0.5```

This will give you two files output_alpha_0p5_path.txt (the path) and  output_alpha_0p5_deg_cross.txt (contains a list of the degree and the crossing number of each vertex). 



# Lowsum

Same as LowMax but the path the algorithm outputs has a low sum of crossing numbers

