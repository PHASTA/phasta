        block data endata
c
c----------------------------------------------------------------------
c 
c  Almost all data statements are stored in this block data.
c
c
c Farzin Shakib, Summer 1985.
c Zdenek Johan,  Winter 1991.  (Fortran 90)
c----------------------------------------------------------------------
c
        include "common.h"
c
c----------------------------------------------------------------------
c
c.... quadrature point data
c
c----------------------------------------------------------------------

	data master / 0 /
c
c.... boundary nodes of boundary elements
c.... common /bndnod/ mnodeb(9,8,3)
c
        data (((mnodeb(i,j,k), i=1,9), j=1,8), k=1,2)
     &              /  1,  0,  0,   0,  0,  0,   0,  0,  0, 
     &                 1,  0,  0,   0,  0,  0,   0,  0,  0, 
     &                 1,  0,  0,   0,  0,  0,   0,  0,  0, 
     &                 1,  0,  0,   0,  0,  0,   0,  0,  0, 
     &                 1,  0,  0,   0,  0,  0,   0,  0,  0, 
     &                 1,  0,  0,   0,  0,  0,   0,  0,  0, 
     &                 1,  0,  0,   0,  0,  0,   0,  0,  0, 
     &                 1,  0,  0,   0,  0,  0,   0,  0,  0,     ! 1D
     &                 1,  2,  0,   0,  0,  0,   0,  0,  0, 
     &                 1,  2,  0,   0,  0,  0,   0,  0,  0, 
     &                 1,  2,  0,   0,  0,  0,   0,  0,  0, 
     &                 1,  2,  0,   0,  0,  0,   0,  0,  0, 
     &                 1,  2,  5,   0,  0,  0,   0,  0,  0, 
     &                 1,  2,  4,   0,  0,  0,   0,  0,  0,
     &                 1,  2,  4,   0,  0,  0,   0,  0,  0,
     &                 1,  2,  4,   0,  0,  0,   0,  0,  0   /  ! 2D
        data (((mnodeb(i,j,k), i=1,9), j=1,8), k=3,3)
     &              /  1,  2,  3,   4,  0,  0,   0,  0,  0, 
     &                 1,  2,  3,   0,  0,  0,   0,  0,  0, 
     &                 1,  2,  3,   0,  0,  0,   0,  0,  0, 
     &                 1,  2,  5,   4,  0,  0,   0,  0,  0, 
     &                 1,  2,  3,   4,  9, 10,  11, 12, 21,
     &                 1,  2,  3,   5,  6,  9,   0,  0,  0,
     &                 1,  2,  3,   7,  9,  8,   0,  0,  0,
     &                 1,  2,  5,   4,  7, 10,  13, 14, 16   /  ! 3D  
c
c----------------------------------------------------------------------
c
c.... element category information
c.... common /elmcat/ mcsyst, melCat, nenCat(8,3), nfaCat(8,3)
c
        data mcsyst, melCat
     &     /   4,      8    /           ! caution: see below
c
        data nenCat /  2,  2,  2,  2,    3,  3,  3,  3,         ! 1D
     &                 4,  3,  3,  4,    9,  6,  6,  9,         ! 2D
     &                 8,  4,  6,  6,   27, 10, 18, 18    /     ! 3D
c
        data nfaCat /  2,  2,  2,  2,    2,  2,  2,  2,         ! 1D
     &                 4,  3,  3,  4,    4,  3,  3,  4,         ! 2D
     &                 6,  4,  5,  5,    6,  4,  5,  5    /     ! 3D
c
c melCat affects: nenCat, nfaCat, mnodeb
c
c----------------------------------------------------------------------
c
c.... maximum number of quadrature points per nsd
c.... common /intpar/ intmax
c
        data intmax /  3  /
c
c----------------------------------------------------------------------
c
c.... io channels
c.... common /io    / iin,    igeom,  ipar,   ibndc,  imat,   iecho,
c....                 iout,   ichmou, irstin, irstou, ihist,  iflux,
c....                 ierror, itable, iforce, igraph, itime
c
        data    iin,    igeom,  ipar,   ibndc,  imat,   iecho,
     &          iout,   ichmou, irstin, irstou, ihist,  iflux,
     &          ierror, itable, iforce, igraph, itime
     &  /       10,     11,     12,     13,     14,     15,
     &          16,     17,     18,     19,     20,     21,
     &          22,     23,     24,     25,     26      /
c
c----------------------------------------------------------------------
c
c.... io file names
c.... common /ioname/ fin,    fgeom,  fpar,   fbndc,  fmat,   fecho,
c....                 frstin, frstou, fhist,  ferror, ftable, fforce,
c....                 ftime
c
        data    fin,    fgeom,  fpar,   fbndc,  fmat,   fecho,
     &          frstin, frstou, fhist,  ferror, ftable, fforce, fgraph,
     &          ftime
     &  /       'input.dat',            'geombc.dat',
     &          'partition.dat',        'bc.dat',
     &          'material.dat',         'echo.dat',
     &          'restart',           'restart',
     &          'histor.dat',           'error.dat',
     &          'table.dat',            'forces.dat',
     &          'graph.dat',            'time.out'       /
c
c----------------------------------------------------------------------
c
c.... run parameters 
c.... common /matpar/ ithm,   pr,     Planck, Stefan, Nh,     Rh,
c                     Rgas,   gamma,  gamma1, s0,     const,  xN2,
c                     xO2,    yN2,    yO2,    Msh,    cpsh,   s0sh,
c                     h0sh,   Rs,     cps,    cvs,    h0s,    Trot,
c                     sigs,   Tvib,   g0s,    dofs
c
c        data    pr
c     &  /       7.20000000000000d-1      /
c
        data    Planck,               Stefan,
     &          Nh,                   Rh,
     &          gamma,                gamma1
     &  /       6.62617600000000d-34, 5.66970000000000d-08,
     &          6.02204500000000d+23, 8.31441000000000d+0,
     &          1.40000000000000d+0,  0.40000000000000d+0    /
c
        data    xN2,                  xO2
     &  /       0.79000000000000d+0,  0.21000000000000d+0     /
c
        data    Msh 
     &  /       2.80000000000000d-2,  3.20000000000000d-2,
     &          3.00000000000000d-2,  1.40000000000000d-2,
     &          1.60000000000000d-2    /
c
        data    h0sh
     &  /       0.00000000000000d+0,  0.00000000000000d+0,
     &          8.97750000000000d+4,  4.70820000000000d+5,
     &          2.46790000000000d+5    /
c
        data    Trot
     &  /       2.87000000000000d+0,  2.08000000000000d+0,
     &          2.45000000000000d+0,  0.00000000000000d+0,
     &          0.00000000000000d+0    /
c
        data    sigs
     &  /       2.00000000000000d+0,  2.00000000000000d+0,
     &          1.00000000000000d+0,  0.00000000000000d+0,
     &          0.00000000000000d+0    /
c
        data    Tvib
     &  /       3.39350000000000d+3,  2.27356000000000d+3,
     &          2.73887000000000d+3,  0.00000000000000d+0,
     &          0.00000000000000d+0    /
c
        data    g0s
     &  /       1.00000000000000d+0,  3.00000000000000d+0,
     &          4.00000000000000d+0,  4.00000000000000d+0,
     &          9.00000000000000d+0    /
c
        data    dofs
     &  /       5.00000000000000d+0,  5.00000000000000d+0,
     &          5.00000000000000d+0,  3.00000000000000d+0,
     &          3.00000000000000d+0    /
c
c----------------------------------------------------------------------
c
c.... dynamic storage pointer management data
c.... common /point / mbeg,   mend,   mprec
c
        data    mbeg,    mend,    mprec
     &  /       1,      100000,     2   /
c
c----------------------------------------------------------------------
c
c.... residual statistics data
c.... common /resdat/ resfrt
c
        data    resfrt
     &  /       0.00000000000000d+0     /
c
c----------------------------------------------------------------------
c
c.... symmetric storage pointers
c.... common /sympar/ indsym(5,5)
c
        data indsym /    1,  2,  4,  7, 11,
     &                   2,  3,  5,  8, 12,
     &                   4,  5,  6,  9, 13,
     &                   7,  8,  9, 10, 14,
     &                  11, 12, 13, 14, 15   /
c
c----------------------------------------------------------------------
c
c.... timer parameters
c.... common /timer1/ ccode(13)
c.... common /timer2/ icd
c
        data    ccode
     &  /       'Input   ', 'PrProces', 'Rezoning', 'Elm_Form',
     &          'Solver  ', 'Bnd_Flux', 'Output  ', 'Mapping ',
     &          'Gather  ', 'Scatter ', 'Begin   ', 'End     ',
     &          'Back    ' /
c
        data    icd
     &  /       11         /
c
c----------------------------------------------------------------------
c
c.... end
c
        end
