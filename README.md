# Documentation

## Licensing Information

Skyler is licensed under an MIT license <https://opensource.org/licenses/MIT>.

You should have received a copy of the MIT Public License along with Skyler.  If not, please see <https://opensource.org/licenses/MIT>.


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
<yossi.bokor@anu.edu.au> \
<yossi.bokor@sydney.edu.au> \
<yossi@gtfo.top> \
[Personal Webpage](http://gtfo.top/yossi.html) \
and \
Christopher Williams\
<christopher.williams@anu.edu.au>\
## Installation
To install Skyler, run the following in `Julia`:
```julia
using Pkg
Pkg.add("Skyler")
```

## Functionality
- Skyler identifies the coarsest abstract graph structure underlying a point cloud, and modells it. Currently, we are restricted to graphs with linear edges which satisfy conditions detailed below. 


### Obtaining Abstract Structure

As input Skyler accepts a point cloud, as an $d \times n$ array, a parameter $\varepsilon$, an inner radius, an outer radius, and an angle condition. The parameter $\varepsilon$ is the same paramter such that the point cloud is an $\varepsilon$-sample of an embedded graph $|G|$. The radii and angle condition are realted to assumptions on the embedding of $|G|$, and their derivations can be found in (DGC ARTICLE I HOPE). To obtain the abstract graph $G$, run 

```julia
structure = skyler(points, sample_epsilon, angle_condition, inner_radius, outer_radius, vertex_threshold, edge_threshold, out="struct")
```

which returns a data frame `partition` and an array `boundaries`.

The columns of `partition` are as follows: the first is an (arbitrary) index of the strata a point as been assocaited to, the second is the dimension of said strata, and the rest are the coordinates of the point. The array `boundaries` is a list of simplicies given by their vertices (i.e for a vertex it just lists the vertex number, and for edges it is a pair).

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
structure = skyler(points, 0.01, 0.08, 0.1, 1/4, 0.03, 0.09, out="test")
```
which should result in output similar to the following:

```julia
(450×4 DataFrames.DataFrame
│ Row │ x1      │ x2      │ x3          │ x4           │
│     │ Float64 │ Float64 │ Float64     │ Float64      │
├─────┼─────────┼─────────┼─────────────┼──────────────┤
│ 1   │ 1.0     │ 0.0     │ -0.00956676 │ -0.00280761  │
│ 2   │ 1.0     │ 0.0     │ 0.00833054  │ -0.000256456 │
│ 3   │ 1.0     │ 0.0     │ 0.0151465   │ -9.84807e-5  │
│ 4   │ 1.0     │ 0.0     │ 0.0105458   │ 0.0150875    │
│ 5   │ 1.0     │ 0.0     │ 0.02135     │ 0.00899252   │
│ 6   │ 1.0     │ 0.0     │ 0.0246492   │ 0.0116247    │
│ 7   │ 1.0     │ 0.0     │ 0.0262769   │ 0.0153962    │
│ 8   │ 1.0     │ 0.0     │ 0.0332015   │ 0.0117242    │
│ 9   │ 1.0     │ 0.0     │ 0.0378463   │ 0.0242734    │
│ 10  │ 1.0     │ 0.0     │ 0.0392313   │ 0.0289288    │
│ 11  │ 1.0     │ 0.0     │ 0.0501769   │ 0.0161012    │
│ 12  │ 1.0     │ 0.0     │ 0.0515084   │ 0.0200041    │
│ 13  │ 1.0     │ 0.0     │ 0.0562974   │ 0.0322489    │
│ 14  │ 1.0     │ 0.0     │ 0.0574664   │ 0.0355348    │
│ 15  │ 1.0     │ 0.0     │ 0.0643197   │ 0.034779     │
│ 16  │ 1.0     │ 0.0     │ 0.0673927   │ 0.0370305    │
│ 17  │ 1.0     │ 0.0     │ 0.0743087   │ 0.0348942    │
│ 18  │ 1.0     │ 0.0     │ 0.0804208   │ 0.0314756    │
│ 19  │ 2.0     │ 0.0     │ 1.9932      │ 0.997527     │
│ 20  │ 2.0     │ 0.0     │ 1.96634     │ 0.993213     │
│ 21  │ 2.0     │ 0.0     │ 1.9432      │ 0.979328     │
│ 22  │ 2.0     │ 0.0     │ 1.92279     │ 0.971437     │
│ 23  │ 2.0     │ 0.0     │ 1.9147      │ 0.954751     │
│ 24  │ 2.0     │ 0.0     │ 1.91975     │ 0.956415     │
⋮
│ 426 │ 3.0     │ 1.0     │ 1.80267     │ 0.90083      │
│ 427 │ 3.0     │ 1.0     │ 1.80532     │ 0.905353     │
│ 428 │ 3.0     │ 1.0     │ 1.8063      │ 0.909872     │
│ 429 │ 3.0     │ 1.0     │ 1.81288     │ 0.912725     │
│ 430 │ 3.0     │ 1.0     │ 1.8159      │ 0.914215     │
│ 431 │ 3.0     │ 1.0     │ 1.82205     │ 0.920487     │
│ 432 │ 3.0     │ 1.0     │ 1.82698     │ 0.917841     │
│ 433 │ 3.0     │ 1.0     │ 1.83481     │ 0.908489     │
│ 434 │ 3.0     │ 1.0     │ 1.84158     │ 0.913923     │
│ 435 │ 3.0     │ 1.0     │ 1.83749     │ 0.925785     │
│ 436 │ 3.0     │ 1.0     │ 1.84788     │ 0.921541     │
│ 437 │ 3.0     │ 1.0     │ 1.85204     │ 0.921254     │
│ 438 │ 3.0     │ 1.0     │ 1.85401     │ 0.930794     │
│ 439 │ 3.0     │ 1.0     │ 1.85932     │ 0.92983      │
│ 440 │ 3.0     │ 1.0     │ 1.86424     │ 0.93014      │
│ 441 │ 3.0     │ 1.0     │ 1.87317     │ 0.926241     │
│ 442 │ 3.0     │ 1.0     │ 1.87273     │ 0.939438     │
│ 443 │ 3.0     │ 1.0     │ 1.87598     │ 0.937406     │
│ 444 │ 3.0     │ 1.0     │ 1.88738     │ 0.9346       │
│ 445 │ 3.0     │ 1.0     │ 1.88842     │ 0.9393       │
│ 446 │ 3.0     │ 1.0     │ 1.89089     │ 0.949122     │
│ 447 │ 3.0     │ 1.0     │ 1.89864     │ 0.94073      │
│ 448 │ 3.0     │ 1.0     │ 1.90089     │ 0.951265     │
│ 449 │ 3.0     │ 1.0     │ 1.90222     │ 0.957706     │
│ 450 │ 3.0     │ 1.0     │ 1.90771     │ 0.95995      │, Any[1, 2, Any[1, 2]], (Array{Any,1}[[1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19], [2, 443, 437, 433, 430, 431, 432, 434, 435, 436  …  440, 441, 442, 444, 445, 446, 447, 448, 449, 450]], 2, Array{Any,1}[[20, 21, 22, 23, 24, 25, 26, 27, 28, 29  …  420, 421, 422, 423, 424, 425, 426, 427, 428, 429]], 1))
```

