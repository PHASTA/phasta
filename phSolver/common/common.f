      subroutine initphcommonvars() bind(C, name='initPhastaCommonVars')
      use iso_c_binding
      include "common.h"

      character(8), dimension(13) :: names
     &   = (/ 'Input   ', 'PrProces', 'Rezoning', 'Elm_Form',
     &        'Solver  ', 'Bnd_Flux', 'Output  ', 'Mapping ',
     &        'Gather  ', 'Scatter ', 'Begin   ', 'End     ',
     &        'Back    ' /)
      ccode = reshape(names,shape(names))

      intmax = 3
      master = 0
      icd = 11

      indsym = reshape((/ 1,  2,  4,  7, 11,
     &                    2,  3,  5,  8, 12,
     &                    4,  5,  6,  9, 13,
     &                    7,  8,  9, 10, 14,
     &                    11, 12, 13, 14, 15   /),
     &                    shape(indsym))

      resfrt = 0.00000000000000d+0

      mbeg = 1
      mend = 100000
      mprec = 2

      fin = 'input.dat'
      fgeom = 'geombc.dat'
      fpar = 'partition.dat'
      fbndc = 'bc.dat'
      fmat = 'material.dat'
      fecho = 'echo.dat'
      frstin = 'restart'
      frstou = 'restart'
      fhist = 'histor.dat'
      ferror = 'error.dat'
      ftable = 'table.dat'
      fforce = 'forces.dat'
      fgraph = 'graph.dat'
      ftime = 'time.out'

      iin = 10
      igeom = 11
      ipar = 12
      ibndc = 13
      imat = 14
      iecho = 15
      iout = 16
      ichmou = 17
      irstin = 18
      irstou = 19
      ihist = 20
      iflux = 21
      ierror = 22
      itable = 23
      iforce = 24
      igraph = 25
      itime = 26

      mcsyst = 4
      melCat = 8
      nenCat = reshape((/  2,  2,  2,  2,    3,  3,  3,  3,        ! 1D
     &                     4,  3,  3,  4,    9,  6,  6,  9,        ! 2D
     &                     8,  4,  6,  6,   27, 10, 18, 18    /),
     &                 shape(nenCat))     ! 3D
      nfaCat = reshape((/  2,  2,  2,  2,    2,  2,  2,  2,       ! 1D
     &                     4,  3,  3,  4,    4,  3,  3,  4,       ! 2D
     &                     6,  4,  5,  5,    6,  4,  5,  5    /), ! 3D
     &                 shape(nfaCat))


      mnodeb = reshape((/  1,  0,  0,   0,  0,  0,   0,  0,  0,
     &                     1,  0,  0,   0,  0,  0,   0,  0,  0,
     &                     1,  0,  0,   0,  0,  0,   0,  0,  0,
     &                     1,  0,  0,   0,  0,  0,   0,  0,  0,
     &                     1,  0,  0,   0,  0,  0,   0,  0,  0,
     &                     1,  0,  0,   0,  0,  0,   0,  0,  0,
     &                     1,  0,  0,   0,  0,  0,   0,  0,  0,
     &                     1,  0,  0,   0,  0,  0,   0,  0,  0,   ! 1D
     &                     1,  2,  0,   0,  0,  0,   0,  0,  0,
     &                     1,  2,  0,   0,  0,  0,   0,  0,  0,
     &                     1,  2,  0,   0,  0,  0,   0,  0,  0,
     &                     1,  2,  0,   0,  0,  0,   0,  0,  0,
     &                     1,  2,  5,   0,  0,  0,   0,  0,  0,
     &                     1,  2,  4,   0,  0,  0,   0,  0,  0,
     &                     1,  2,  4,   0,  0,  0,   0,  0,  0,
     &                     1,  2,  4,   0,  0,  0,   0,  0,  0,  ! 2D
     &                     1,  2,  3,   4,  0,  0,   0,  0,  0,
     &                     1,  2,  3,   0,  0,  0,   0,  0,  0,
     &                     1,  2,  3,   0,  0,  0,   0,  0,  0,
     &                     1,  2,  5,   4,  0,  0,   0,  0,  0,
     &                     1,  2,  3,   4,  9, 10,  11, 12, 21,
     &                     1,  2,  3,   5,  6,  9,   0,  0,  0,
     &                     1,  2,  3,   7,  9,  8,   0,  0,  0,
     &                     1,  2,  5,   4,  7, 10,  13, 14, 16   /),
     &                     shape(mnodeb))  ! 3D


      Planck = 6.62617600000000d-34
      Stefan = 5.66970000000000d-08
      Nh = 6.02204500000000d+23
      Rh = 8.31441000000000d+0
      gamma = 1.40000000000000d+0
      gamma1 = 0.40000000000000d+0
      xN2 = 0.79000000000000d+0
      xO2 = 0.21000000000000d+0
      Msh = reshape((/ 2.80000000000000d-2,  3.20000000000000d-2,
     &                 3.00000000000000d-2,  1.40000000000000d-2,
     &                 1.60000000000000d-2 /),
     &                 shape(Msh))
      h0sh = reshape((/ 0.00000000000000d+0,  0.00000000000000d+0,
     &                  8.97750000000000d+4,  4.70820000000000d+5,
     &                  2.46790000000000d+5 /),
     &                  shape(h0sh))
      Trot = reshape((/ 2.87000000000000d+0,  2.08000000000000d+0,
     &                  2.45000000000000d+0,  0.00000000000000d+0,
     &                  0.00000000000000d+0 /),
     &                  shape(Trot))
      sigs = reshape((/ 2.00000000000000d+0,  2.00000000000000d+0,
     &                  1.00000000000000d+0,  0.00000000000000d+0,
     &                  0.00000000000000d+0 /),
     &                  shape(sigs))
      Tvib = reshape((/ 3.39350000000000d+3,  2.27356000000000d+3,
     &                  2.73887000000000d+3,  0.00000000000000d+0,
     &                  0.00000000000000d+0 /),
     &                  shape(Tvib))
      g0s = reshape((/ 1.00000000000000d+0,  3.00000000000000d+0,
     &                 4.00000000000000d+0,  4.00000000000000d+0,
     &                 9.00000000000000d+0 /),
     &                 shape(g0s))
      dofs = reshape((/ 5.00000000000000d+0,  5.00000000000000d+0,
     &                  5.00000000000000d+0,  3.00000000000000d+0,
     &                  3.00000000000000d+0 /),
     &                  shape(dofs))
      end subroutine initphcommonvars
