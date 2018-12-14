This Finite Element source code is dated 1991-1993 . ie my early post-doc days.
It was my reqworking of my Ph.D software to make it evelove into a more general purpose framework.
ie. moving from just analysing SLope Stability problems to being used to analyse Spud-cans and Piles and later eigenmodes.

Daniel Kidger . daniel.kidger@gmail.com

## Notes
  * The Fortran source, although generic has some embedded compiler hints for the Salford Fortran Compiler for MS-DOS.
  * The original filenames were all upercase with the .FOR extension. These renamed here in a lowercase and with .f extensions
  * This is Fortan77 not Fortran90 so is fixed format. It does though use every more modern Fortran features that compilers alow, eg lowercase names, and no lomit of 6 characters.


## Important Files
| filename   | functions | Description |
| ---------- | --------- | ----------- |
| **Libraries** | 
| funs.f     |   25     | Shape functions and derivatives of all known elements - call generically via shape() 
| shape.f    |   25     | varient?
| l1a.f      |   34     | : main library of matrix algebra and FEA subroutines
|**Programs**|    
| gen4.f     |   8  | Main program for SLope Stability problems  - I think this also calculates bending moments for Piles and Spud-cans
| gen6.f | x | Later vbersion of gen4.f - see revsion history in the source.
| eigen_4.f  |   1  | Creates an inididual Finite Element Stiffness matrix then solves the EIgenvalue problem to show fundamental deformation modes. Compile and link with l1a.f
| pl8d.f |49 | Graphical post processor - forerunner of DANPLOT. It #include's shape.f  Note the built-in PostScript driver :-) .  This is monolythinc with 48 subroutines!

## Data Files

| filename   | run with | Description |
| ---------- | -------- | ----------- |
| spud3d4.d | gen4 ? | Spud can mesh as a quarter cylinder of soil with two piles close to the spud can.
| tun_12.d | gen6 ~? |  Showing use of *DANBLOCKS to create a lined tunnel with a block of ground 

## Other Files
| filename   | functions | Description |
| ---------- | --------- | ----------- |
Gen5d.lst | | Output from the Salford Fortran Compiler /Linker
mesh4.f | | Mesh generator - Uses the *Keyword Syntax just create a mesh but not analyse it. This means we need much less memory footprint as no KV array etc. Note that this has dependances on Salford Fortran's command line function CMNAM@()
mesh4af | |
mesh4b.f | |
mesh4c.f | |
mohr.f1 | 7 | Just a cut out of the Mohr-Coulomb Subroutines
mung.f | 1 | source to source tool to allow use of the Waterloo Fortran compiler (It doesn't like !comments)


