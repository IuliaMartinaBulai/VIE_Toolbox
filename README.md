# VIE_Toolbox
Volterra Integral Equation Toolbox

(c) Developed by Iulia Martina Bulai, Maria Carmela De Bonis, Concetta Laurita, Valeria Sagaria

This toolbox contains code implementing the Volterra Integral Equation solving method, a general framework 
for the numerical resolution of Volterra integral equations of the second kind. 
In this toolbox we propose two main algorithms, (i) one for approximating the solutions of the more
general VIE of second kind and (ii) a second one for computing the biological observables, such as the total 
metastatic mass and the cumulative number of metastases, respectively, from the generalized metastatic tumor
growth model that describes the primary tumor growth by means of an Ordinary Differential Equation (ODE) and 
the evolution of the metastatic density using a transport Partial Differential Equation (PDE), [2]. 

The software implements a numerical method developed in the following two very recent papers: <br>
[1] M.C. De Bonis, C. Laurita, V. Sagaria. A numerical method for linear Volterra integral equations on 
infinite intervals and its application to the resolution of metastatic tumor growth. Appl. Num. Math. 
172 (2022), 475-496. <br>
[2] I.M. Bulai, M.C. De Bonis, C. Laurita, V. Sagaria. Modeling metastatic tumor evolution, 
numerical resolution and growth prediction. submitted 2022.

The algorithms implemented here are described in detail in: 

[3] I.M. Bulai, M.C. De Bonis, C. Laurita, V. Sagaria. MatLab Toolbox for the numerical solution 
of linear Volterra integral equations arising in metastatic tumor growth models. submitted 2022.


GETTING STARTED
Run

&gt;&gt; vie_demo1 (Figure 2 from [3])

&gt;&gt; vie_demo2 (Figure 3 from [3])

&gt;&gt; vie_demo3 (Figure 4 from [3])

&gt;&gt; vie_demo4 (Figure 1 from [2])

&gt;&gt; vie_demo5 (Figure 2 from [2])

as a demonstration of what this toolbox is used for.

-The files argselectAssign.m argselectCheck.m are part of Spectral Graph Wavelet Transform (SGWT) 
toolbox and can be downloaded at page: https://wiki.epfl.ch/sgwt
-The file gaussq.m is the result of a translation of the fortran procedure http://www.netlib.org/go/gaussq.f

License : 

The VIE toolbox is a Matlab library released under the GPL.

The VIE toolbox is free software: you can redistribute it and/or modify it under the terms of the GNU 
General Public License as published by the Free Software Foundation, either version 3 of the License, 
or (at your option) any later version.

The VIE toolbox is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the VIE toolbox.  
If not, see <http://www.gnu.org/licenses/>.
