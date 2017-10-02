---
{title: swplot.arrow, link: swplot.arrow, summary: draws a 3D arrow using patch, keywords: sample,
  sidebar: sw_sidebar, permalink: swplot_arrow, folder: swplot, mathjax: 'true'}

---

### Syntax

`hpatch = swplot.arrow(rstart, rend, r, α, lhead, {npatch})`

### Description

hPatch = SWPLOT.ARROW(handle,...)
 
Handle can be the handle of an axes object or a patch object. It either
selects an axis to plot or a patch object (triangulated) to add vertices
and faces.
 

### Input Arguments

`handle`
: Handle of an axis or patch object. In case of patch object, the
  constructed faces will be added to the existing object instead
  of creating a new one.

`rStart`
: Coordinate of the starting point.

`rEnd`
: Coordinate of the end point.

`R`
: Radius of the arrow body.

`α`
:   Angle of the head in °s.

`lHead`
: Length of the head.

`nPatch`
: Number of points on the curve, default value is stored in
  swpref.getpref('npatch').

### See Also

[swplot.cylinder](swplot_cylinder)

{% include links.html %}
