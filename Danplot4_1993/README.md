# DANPLOT - A 3D Graphical FIniter Element Post Proessor.


This is an instance of Danplot as of c. Nov.1994


## Important Files
| filename   | functions | Description |
| ---------- | --------- | ----------- |
| bits.f | 14 | a sunbset of FE5lib plus some routines that dint firt anywhere else - eg the Danplot Splash Screen 
| danpl.f | x | The main program . aha Plotter_7C
| draw.f 29 | main graphics routines: Gouraud shader, 3D camera positioning, .. Note the PCX bitmap reading routines.
| elemets.f | x | put/get elements intio teh internal data structures
| facet.f | 5 | sub-faceting fstrip2(), get_fac, etc.
| files.f | 3 | data file parser, copes with comments, #includes etc.
| k_mesh.f | 10 | Keyword handler for importing FE mesh file formats.
| k_mmung.f | 26 | DANBLOCKS p;lus a whiole host of other mesh transformation routines.
| makefile | x | For danplot.exe using Salford Fortran compiler
| menu.f | 24 | Post all 18 menus in the GUI
| read.f | 6 | a few routines for reading meshes
| shape.f | 30 | Full set of shape functiuons of all known elements - invcluding the 14node brick by the author :-)
| _driver.f | 0 | test progra, to simply draw the gui menus
| 



## Other Files
| filename   | functions | Description |
| ---------- | --------- | ----------- |
| eigen_4.f | x | Driver program for finite element Stiffness Matrix Eigenmode analysis.
| reame.f | x |  Nop-FEA program BY DR. Y. H. HUANG, PROFESSOR OF CIVIL ENGINEERING, UNIVERSITY OF KENTUCKY. Calculates stability of soil slopes by means of trial slip circles.
| inreame.dat | x | Some data file - not sure what it is.
| tblock.f | 2 | Random nukmber tester? Maybe by Carolus S.?  Plots them using Salford Fortrean hraphics routines.
| tirgeo2.f | x |  _Thameside  Industrial Route embankment_ . Carolus?
