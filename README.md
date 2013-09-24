MOLPOP-CEP
=============

MOLPOP-CEP is a code for the exact solution of radiative transfer problems in multi-level
atomic systems. The novel contribution of the code is that the radiative transfer equations
is analytically integrated so that the final problem is reduced to the solution of a non-
linear algebraic system of equations in the level populations. The radiative transfer is solved
analytically using the Coupled Escape Probability formalism presented by Elitzur & Asensio
Ramos (2006). The current version of
the code is limited to plane-parallel slabs that can present arbitrary spatial variations of the
physical conditions.

The code is written in standard Fortran 90. It is based on the MOLPOP code written by
M. Elitzur that used single-zone escape probabilities for the solution of the radiative transfer
problem. The original MOLPOP code, written in Fortran 77, has been ported to Fortran
90. During the translation, the code has been modularized and all common blocks have been
moved to external modules that are used only where necessary. All the machinery present in
MOLPOP for reading the input file and carry out all the needed calculations (interpolation
of the collisional rates, selection of the active levels, etc.) are still present. The fundamental
idea when merging together the MOLPOP code and the CEP code was to maintain the large
flexibility of the input already present in MOLPOP. When the solution method chosen in
the input file is the single-zone escape probability, the original MOLPOP code is executed.
When CEP is chosen as the method, the routines belonging to the CEP code are used.
Although the resulting code is a mixture of two existing codes, the interface between both
is simple and robust.

Some parts of the code can be run online on 
http://www.iac.es/proyecto/inversion/online/molpop_code/molpop.php
and we are working towards implementing more and more options on the
online version.