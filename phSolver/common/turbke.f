c-----------------------------------------------------------------------
c
c     Lam-Bremhorst Kay-Epsilon turbulence model constants 
c
c-----------------------------------------------------------------------
      module turbKE

      real*8 ke_C_mu, ke_C_eps1, ke_C_eps2
      real*8 ke_sigma, ke_otherstuff
      parameter ( 
     &     ke_C_mu       = 0.09,
     &     ke_C_eps1     = 1.44,
     &     ke_C_eps2     = 1.92,
     &     ke_sigma      = 1.3,
     &     ke_otherstuff = 1.50d0
     &     )

      
      end module

