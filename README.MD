Maximizing point count filter script
====================================

Issue addressed by the script
-----------------------------

Given set of points distributed in geodetic spherical space (EPSG:4326), script produces subset of maximal cardinality (maximal number of points) such that no two or more points are closer than given distance constraint (in meters).


The script takes array of points, stored in structured csv format with following field order: id, longitude, latitude; separated by comma (,) and having field names in the first row, and produces a filtered list in the csv file of the same structure as input file.


Sample set of randomly generated points can be found in the project's root folder.


Dependencies
------------

python ver > 3<br/>
numpy  ver > 1.11.0


Usage
-----

Example:
> python filter.py inputDS outputDS constraining_distance (m)

where:
> inputDS  : filename of input csv file<br>
> outputDS : filename of output csv file<br>
> constraining_distance : minimal distance between any two or more points


Disclaimer
----------

This script has not beet throughly tested nor correctest of the output is guaranteed by the author. The script is possible implementation of solution to given problem and is shipped "AS IS" without warranty of any kind or ant liability from the side of the script's author. The usage and distribution of the script, including results produced by the same, is free without any contribution to the author.