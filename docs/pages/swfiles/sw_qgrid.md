---
{title: sw_qgrid, link: sw_qgrid, summary: creates a Q grid, keywords: sample, sidebar: sw_sidebar,
  permalink: sw_qgrid, folder: swfiles, mathjax: 'true'}

---
  
### Syntax
  
`qgrid = sw_qgrid(Name,Value)`
  
### Description
  
`qgrid = sw_qgrid(Name,Value)` generates n-dimensional grids ($$n<=3$$) in
3D space, e.g. points on a line in 3D. It uses $$n$$ linearly independent
vectors ("lattice vectors") and bin values (coordinates in "lattice
units" or "lu") to generate the points. It works similarly as the d3d
constructor in [Horace](http://horace.isis.rl.ac.uk/Main_Page).
  
### Name-Value Pair Arguments
  
`'u'`
:  Row vector with 3 elements, determines the first axis in 3D
   space, default value is `[1 0 0]`.
  
`'v'`
:  Second axis, default value is `[0 1 0]`.
  
`'w'`
:  Third axis, default value is `[0 0 1]`.
  
`'uoffset'`
:  Row vector with 3 elements, determines the offset of origin
   in lu, (fourth element is accepted but discarded).
  
`'ubin'`
:  Bin points along the first axis. Can be a vector with 1, 2 or 3
   elements:
 
   * `[B1]`        single value along the $$u$$-axis at a coordinate of `B1*u`
   * `[B1 B2]`     range along the $$u$$-axis at coordinates of `[B1:1/nExt:B2]*u`
   * `[B1 dB B2]`  range along the $$u$$-axis at coordinates of `[B1:dB:B2]*u`
  
`'vbin'`
:  Same as `ubin` but along the $$v$$-axis.
  
`'wbin'`
:  Same as `ubin` but along the $$w$$-axis.
  
`'nExt'`
:  Vector with $$n$$-elements that can define fractional bin steps,
   default values is `[1 1 1]`.
  
`'lab'`
:  Cell array of projection axis labels with 3 elements (4th
   element discarded), e.g. `{'x' 'y' 'z'}`.
  
The dimension count $$n$$ is determined by the number of given bins
($$1<=n<=3$$), so if only `ubin` is given, $$n=1$$; if both `ubin` and `vbin`
are defined then $$n=2$$, etc.
  
### Output Arguments
  
`qGrid`
: A matrix with dimensions of $$[3\times n_{ax1}\times n_{ax2},...]$$,
  where $$n_{axi}$$ is the index of points along $$i$$th axis with $$1<=i<=n$$.
  
### See Also
  
[sw_qscan](sw_qscan)
 

{% include links.html %}
