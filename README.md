# Documentation

## Licensing Information

Skyler is licensed under an MIT license <https://opensource.org/licenses/MIT>.

You should have received a copy of the MIT Public License along with Skyler.  If not, please see <https://opensource.org/licenses/MIT>.
[![DOI](https://zenodo.org/badge/252151758.svg)](https://zenodo.org/badge/latestdoi/252151758)


```
Begin license text.
Skyler
Copyright (C) 2020 Yossi Bokor

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

End license text.
```

Skyler is free software: you can redistribute it and/or modify
 it under the terms of the MIT License as published by
 the Open Source Initiative

 Skyler is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the MIT License for more details.

## Contributors

Skyler is produced and maintained by

Yossi Bokor \
<yossi@yossi.eu> \
[Personal Webpage](http://yossi.eu) \
and \
Christopher Williams\
<christopher.williams@anu.edu.au>
## Installation
To install Skyler, run the following in `Julia`:
```julia
using Pkg
Pkg.add("Skyler")
```

## Functionality
- Skyler identifies the coarsest abstract graph structure underlying a point cloud, and modells it. Currently, we are restricted to graphs with linear edges which satisfy conditions detailed below. 
- You can read the article [Reconstructing linearly embedded graphs: A first step to stratified space learning](https://www.aimsciences.org/article/doi/10.3934/fods.2021026) which introduces the algorithm used in Skyler.
 </ul>


### Obtaining Abstract Structure

As input Skyler accepts a point cloud, as an $d \times n$ array, a parameter $\varepsilon$, an inner radius, an outer radius, and an angle condition. The parameter $\varepsilon$ is the same paramter such that the point cloud is an $\varepsilon$-sample of an embedded graph $|G|$. The radii and angle condition are realted to assumptions on the embedding of $|G|$, and their derivations can be found in (DGC ARTICLE I HOPE). To obtain the abstract graph $G$, run 

```julia
ppc, model_verts, avg_verts = skyler(points, sample_epsilon, radius, EM_it=3, sigma=sample_epsilon/2)
```

which returns a PartitionedPointCloud  `ppc`, an array `model_verts` containing the optimized vertex locations, and an array `avg_verts` containing the average point of each vertex cluster..

A `PartitionedPointCloud` has the following fields:
	- `Points` which is an array of all the points,
	- `Strata` a `Dict{Int, Set}` which collates which points have been assigned to which stratum
	- `Dimensions` a `Dict{Int, Set}` listing which strata are of each dimesnions, 
	- `Boundaries` an array which represents the boundary operator.
	

### Modeling the Embedding

To model the underlying structure, use 

```julia
vertex_locations = structure= skyler(points, sample_epsilon, angle_condition, inner_radius, outer_radius, vertex_threshold, edge_threshold)
```

which returns an $d \times n_v$ array, where $n_v$ is the number of vertices detected, with each column a modelled vertex location.

## Examples

Skyler comes with some point clouds for you can work through as examples. To load the points for each example, run

```julia
points = Load_Example(i)
```
where `i` is the example number.

#### Example 1

Example 1 is a point cloud sampled from a line segment. Load the sample using 

```julia
points = Load_Example(1)
```

Then run Skyler to obtaint the abstract structure and partition by executing 

```julia
ppc = skyler(points, 0.01, 0.12, out="PPC")
```
which should result in output similar to the following:

```julia
The strata are Dict{Int64,Set}(0 => Set(Any[2, 1]),1 => Set(Any[3])) and the boundary matrix is  [0.0 0.0 1.0; 0.0 0.0 1.0; 0.0 0.0 0.0].
PartitionedPointCloud([-0.009566763308567242 1.9932036826253814 … 1.999301914407944 2.003116429018376; -0.0028076084902411745 0.9975271026710401 … 0.9951311441125332 0.99744890927049], Dict{Int64,Set}(2 => Set(Any[2, 441, 442, 432, 447, 428, 430, 437, 435, 431  …  434, 429, 444, 436, 446, 439, 440, 450, 438, 449]),3 => Set(Any[288, 306, 29, 300, 289, 74, 176, 57, 285, 318  …  341, 186, 321, 420, 423, 271, 23, 315, 322, 218]),1 => Set(Any[12, 4, 18, 3, 16, 11, 5, 21, 20, 7, 9, 13, 10, 14, 19, 17, 8, 15, 6, 1])), Dict{Int64,Set}(0 => Set(Any[2, 1]),1 => Set(Any[3])), [0.0 0.0 1.0; 0.0 0.0 1.0; 0.0 0.0 0.0])
```


#### Example 2 

Example 2 is a point cloud sampled from a graph with 5 edges and 5 vertices. Load the sample using 

```julia
points = Load_Example(2)
```

Then run Skyler to obtaint the abstract structure and partition by executing 

```julia
 ppc = skyler(points, 0.1, 1.2, out="PPC")
```
which should result in output similar to the following:

```julia
The strata are Dict{Int64,Set}(0 => Set(Any[4, 2, 3, 5, 1]),1 => Set(Any[7, 9, 10, 8, 6])) and the boundary matrix is  [0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0; 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0].
PartitionedPointCloud([-0.007670121199078972 -0.09959009503870764 … 4.8232765059435225 4.923038846261083; -0.007670121199078972 -0.09959009503870764 … 0.5335914379931135 0.5214683517413553; -0.007670121199078972 -0.09959009503870764 … 3.4192839435193365 3.4814584770681027], Dict{Int64,Set}(7 => Set(Any[241, 197, 215, 249, 207, 201, 283, 252, 182, 279  …  268, 281, 243, 191, 222, 277, 271, 255, 218, 276]),4 => Set(Any[288, 306, 300, 296, 428, 289, 435, 20, 285, 448  …  24, 429, 427, 446, 439, 23, 305, 438, 449, 301]),9 => Set(Any[532, 520, 491, 478, 542, 499, 477, 509, 494, 521  …  519, 560, 540, 535, 562, 485, 502, 498, 496, 508]),10 => Set(Any[633, 658, 654, 624, 611, 614, 625, 612, 616, 664  …  629, 666, 667, 646, 663, 657, 640, 676, 661, 659]),2 => Set(Any[461, 11, 464, 462, 8, 323, 458, 318, 459, 308  …  319, 456, 321, 454, 312, 317, 463, 472, 315, 322]),3 => Set(Any[584, 574, 698, 699, 566, 582, 573, 569, 694, 14  …  576, 688, 695, 578, 689, 580, 687, 686, 15, 581]),5 => Set(Any[148, 136, 25, 147, 29, 151, 144, 155, 142, 150  …  26, 146, 138, 145, 28, 149, 27, 137, 141, 30]),8 => Set(Any[329, 370, 365, 391, 400, 342, 384, 375, 372, 407  …  415, 341, 378, 389, 420, 423, 424, 358, 349, 405]),6 => Set(Any[89, 134, 131, 74, 57, 78, 112, 70, 106, 121  …  81, 98, 51, 73, 119, 53, 116, 123, 56, 108]),1 => Set(Any[47, 32, 2, 40, 587, 171, 39, 46, 158, 43  …  5, 45, 163, 168, 588, 603, 164, 602, 41, 1])…), Dict{Int64,Set}(0 => Set(Any[4, 2, 3, 5, 1]),1 => Set(Any[7, 9, 10, 8, 6])), [0.0 0.0 … 0.0 1.0; 0.0 0.0 … 1.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0])
```
