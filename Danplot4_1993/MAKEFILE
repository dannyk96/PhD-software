#
#
#  This Makefile will compile Dan Kidger's Finite Element Menu Prog
#      15-6-93
#
# f77debug= /fullcheck /undef

#       .. DANPLOT routines ...
#                          .... and some DANLIB subroutines too
objs=  danpl.obj    \
       menu.obj     \
       read.obj     \
       draw.obj     \
       bits.obj     \
       facet.obj    \
       k_mesh.obj   \
       k_mmung.obj  \
       shape.obj    \
       elements.obj \
       files.obj     
                  
.for.obj:
     ftn77  $<  $(f77debug)

danplot.exe: $(objs)
     link77  makefile.lnk
     copy danplot.exe c:\util
     @echo ++++ 
     @echo ++++  DANPLOT has now been compiled and linked :-)
     @echo ++++ 

clean:
     rm *.obj *.lib

#
#
#

f77debug= /fullcheck /undef

objs= test2.obj

.f.obj: 
     ftn77  $< /intl /dreal $(f77debug)

prog.exe: $(objs)
   link77 makefile.lnk


