      subroutine Spsi3pre (BDiag,  lhsK,  col, row)
c
c------------------------------------------------------------------------------
c This is the initialization routine for the Sparse-GMRES solver.
c It pre-preconditions the LHS mass matrix and sets up the 
c sparse preconditioners. (pre-preconditioning is block diagonal
c scaling). 
c
c input:
c     BDiag  (nshg,nflow,nflow)    : block diagonal scaling matrix 
c                                  which is already LU factored.
c     lhsK(nflow*nflow,nnz_tot) : sparse LHS mass matrix
c
c output:
c     lhsK(nflow*nflow,nnz_tot)  : pre-preconditioned (scaled) mass matrix
c
c Nahid Razmara, Spring 2000.	(Sparse Matrix)
c------------------------------------------------------------------------------
c
      use pointer_data

      include "common.h"
c
        integer col(nshg+1),  row(nnz*nshg)
        integer sparseloc, c, c1
        real*8 lhsK(nflow*nflow,nnz_tot)

c
      
      dimension  BDiag(nshg,nflow,nflow)

      
c
c.... block-diagonal pre-precondition LHS 
c
c     
c.... reduce by columns, (left block diagonal preconditioning)
c
c     lhsK  <-- inverse(L^tilde) lhsK     
c     
c
      if(nflow.ne.5) then
         write(*,*)' spsi3pre.f assumed nflow=5'
         stop
      endif
        do i = 1, nshg 
c
            do k = col(i), col(i+1)-1
c     
                  lhsK(2,k) = lhsK(2,k) 
     &                 - BDiag(i,2,1) * lhsK(1,k)
                  lhsK(7,k) = lhsK(7,k) 
     &                 - BDiag(i,2,1) * lhsK(6,k)
                  lhsK(12,k) = lhsK(12,k) 
     &                 - BDiag(i,2,1) * lhsK(11,k)
                  lhsK(17,k) = lhsK(17,k) 
     &                 - BDiag(i,2,1) * lhsK(16,k)
                  lhsK(22,k) = lhsK(22,k) 
     &                 - BDiag(i,2,1) * lhsK(21,k)
c     
                  lhsK(3,k) = lhsK(3,k) 
     &                 - BDiag(i,3,1) * lhsK(1,k)
     &                 - BDiag(i,3,2) * lhsK(2,k)
                  lhsK(8,k) = lhsK(8,k) 
     &                 - BDiag(i,3,1) * lhsK(6,k)
     &                 - BDiag(i,3,2) * lhsK(7,k)
                  lhsK(13,k) = lhsK(13,k) 
     &                 - BDiag(i,3,1) * lhsK(11,k)
     &                 - BDiag(i,3,2) * lhsK(12,k)
                  lhsK(18,k) = lhsK(18,k) 
     &                 - BDiag(i,3,1) * lhsK(16,k)
     &                 - BDiag(i,3,2) * lhsK(17,k)
                  lhsK(23,k) = lhsK(23,k) 
     &                 - BDiag(i,3,1) * lhsK(21,k)
     &                 - BDiag(i,3,2) * lhsK(22,k)
c     
                  lhsK(4,k) = lhsK(4,k)
     &                 - BDiag(i,4,1) * lhsK(1,k) 
     &                 - BDiag(i,4,2) * lhsK(2,k) 
     &                 - BDiag(i,4,3) * lhsK(3,k)
                  lhsK(9,k) = lhsK(9,k)
     &                 - BDiag(i,4,1) * lhsK(6,k) 
     &                 - BDiag(i,4,2) * lhsK(7,k) 
     &                 - BDiag(i,4,3) * lhsK(8,k)
                  lhsK(14,k) = lhsK(14,k)
     &                 - BDiag(i,4,1) * lhsK(11,k) 
     &                 - BDiag(i,4,2) * lhsK(12,k) 
     &                 - BDiag(i,4,3) * lhsK(13,k)
                  lhsK(19,k) = lhsK(19,k)
     &                 - BDiag(i,4,1) * lhsK(16,k) 
     &                 - BDiag(i,4,2) * lhsK(17,k) 
     &                 - BDiag(i,4,3) * lhsK(18,k)
                  lhsK(24,k) = lhsK(24,k)
     &                 - BDiag(i,4,1) * lhsK(21,k) 
     &                 - BDiag(i,4,2) * lhsK(22,k) 
     &                 - BDiag(i,4,3) * lhsK(23,k)
c     
                  lhsK(5,k) = lhsK(5,k)
     &                 - BDiag(i,5,1) * lhsK(1,k) 
     &                 - BDiag(i,5,2) * lhsK(2,k) 
     &                 - BDiag(i,5,3) * lhsK(3,k) 
     &                 - BDiag(i,5,4) * lhsK(4,k) 
                  lhsK(10,k) = lhsK(10,k)
     &                 - BDiag(i,5,1) * lhsK(6,k) 
     &                 - BDiag(i,5,2) * lhsK(7,k) 
     &                 - BDiag(i,5,3) * lhsK(8,k) 
     &                 - BDiag(i,5,4) * lhsK(9,k) 
                  lhsK(15,k) = lhsK(15,k)
     &                 - BDiag(i,5,1) * lhsK(11,k) 
     &                 - BDiag(i,5,2) * lhsK(12,k) 
     &                 - BDiag(i,5,3) * lhsK(13,k) 
     &                 - BDiag(i,5,4) * lhsK(14,k) 
                  lhsK(20,k) = lhsK(20,k)
     &                 - BDiag(i,5,1) * lhsK(16,k) 
     &                 - BDiag(i,5,2) * lhsK(17,k) 
     &                 - BDiag(i,5,3) * lhsK(18,k) 
     &                 - BDiag(i,5,4) * lhsK(19,k) 
                  lhsK(25,k) = lhsK(25,k)
     &                 - BDiag(i,5,1) * lhsK(21,k) 
     &                 - BDiag(i,5,2) * lhsK(22,k) 
     &                 - BDiag(i,5,3) * lhsK(23,k) 
     &                 - BDiag(i,5,4) * lhsK(24,k) 
               enddo
            enddo
            
        do i = 1, nshg
c
            do k = col(i), col(i+1)-1
                j = row(k)
c     
c.... reduce by rows, (right block diagonal preconditioning)
c
c     lhsK   <-- lhsK  inverse(U^tilde)
c     


                       lhsK(1,k)  = BDiag(j,1,1)*lhsK(1,k) 
                       lhsK(2,k)  = BDiag(j,1,1)*lhsK(2,k) 
                       lhsK(3,k)  = BDiag(j,1,1)*lhsK(3,k) 
                       lhsK(4,k)  = BDiag(j,1,1)*lhsK(4,k) 
                       lhsK(5,k)  = BDiag(j,1,1)*lhsK(5,k) 
      
                  lhsK(6,k) = BDiag(j,2,2)*( lhsK(6,k) 
     &                 - BDiag(j,1,2) * lhsK(1,k) )
                  lhsK(7,k) = BDiag(j,2,2)*( lhsK(7,k) 
     &                 - BDiag(j,1,2) * lhsK(2,k) )
                  lhsK(8,k) = BDiag(j,2,2)*( lhsK(8,k) 
     &                 - BDiag(j,1,2) * lhsK(3,k) )
                  lhsK(9,k) = BDiag(j,2,2)*( lhsK(9,k) 
     &                 - BDiag(j,1,2) * lhsK(4,k) )
                  lhsK(10,k) = BDiag(j,2,2)*( lhsK(10,k) 
     &                 - BDiag(j,1,2) * lhsK(5,k) )
c     
                  lhsK(11,k) = BDiag(j,3,3)*( lhsK(11,k) 
     &                 - BDiag(j,1,3) * lhsK(1,k)
     &                 - BDiag(j,2,3) * lhsK(6,k) )
                  lhsK(12,k) = BDiag(j,3,3)*( lhsK(12,k) 
     &                 - BDiag(j,1,3) * lhsK(2,k)
     &                 - BDiag(j,2,3) * lhsK(7,k) )
                  lhsK(13,k) = BDiag(j,3,3)*( lhsK(13,k) 
     &                 - BDiag(j,1,3) * lhsK(3,k)
     &                 - BDiag(j,2,3) * lhsK(8,k) )
                  lhsK(14,k) = BDiag(j,3,3)*( lhsK(14,k) 
     &                 - BDiag(j,1,3) * lhsK(4,k)
     &                 - BDiag(j,2,3) * lhsK(9,k) )
                  lhsK(15,k) = BDiag(j,3,3)*( lhsK(15,k) 
     &                 - BDiag(j,1,3) * lhsK(5,k)
     &                 - BDiag(j,2,3) * lhsK(10,k) )
c     
                  lhsK(16,k) = BDiag(j,4,4)*( lhsK(16,k)
     &                 - BDiag(j,1,4) * lhsK(1,k) 
     &                 - BDiag(j,2,4) * lhsK(6,k) 
     &                 - BDiag(j,3,4) * lhsK(11,k) )
                  lhsK(17,k) = BDiag(j,4,4)*( lhsK(17,k)
     &                 - BDiag(j,1,4) * lhsK(2,k) 
     &                 - BDiag(j,2,4) * lhsK(7,k) 
     &                 - BDiag(j,3,4) * lhsK(12,k) )
                  lhsK(18,k) = BDiag(j,4,4)*( lhsK(18,k)
     &                 - BDiag(j,1,4) * lhsK(3,k) 
     &                 - BDiag(j,2,4) * lhsK(8,k) 
     &                 - BDiag(j,3,4) * lhsK(13,k) )
                  lhsK(19,k) = BDiag(j,4,4)*( lhsK(19,k)
     &                 - BDiag(j,1,4) * lhsK(4,k) 
     &                 - BDiag(j,2,4) * lhsK(9,k) 
     &                 - BDiag(j,3,4) * lhsK(14,k) )
                  lhsK(20,k) = BDiag(j,4,4)*( lhsK(20,k)
     &                 - BDiag(j,1,4) * lhsK(5,k) 
     &                 - BDiag(j,2,4) * lhsK(10,k) 
     &                 - BDiag(j,3,4) * lhsK(15,k) )
c     
                  lhsK(21,k) = BDiag(j,5,5)*( lhsK(21,k)
     &                 - BDiag(j,1,5) * lhsK(1,k) 
     &                 - BDiag(j,2,5) * lhsK(6,k) 
     &                 - BDiag(j,3,5) * lhsK(11,k) 
     &                 - BDiag(j,4,5) * lhsK(16,k) )
                  lhsK(22,k) = BDiag(j,5,5)*( lhsK(22,k)
     &                 - BDiag(j,1,5) * lhsK(2,k) 
     &                 - BDiag(j,2,5) * lhsK(7,k) 
     &                 - BDiag(j,3,5) * lhsK(12,k) 
     &                 - BDiag(j,4,5) * lhsK(17,k) )
                  lhsK(23,k) = BDiag(j,5,5)*( lhsK(23,k)
     &                 - BDiag(j,1,5) * lhsK(3,k) 
     &                 - BDiag(j,2,5) * lhsK(8,k) 
     &                 - BDiag(j,3,5) * lhsK(13,k) 
     &                 - BDiag(j,4,5) * lhsK(18,k) )
                  lhsK(24,k) = BDiag(j,5,5)*( lhsK(24,k)
     &                 - BDiag(j,1,5) * lhsK(4,k) 
     &                 - BDiag(j,2,5) * lhsK(9,k) 
     &                 - BDiag(j,3,5) * lhsK(14,k) 
     &                 - BDiag(j,4,5) * lhsK(19,k) )
                  lhsK(25,k) = BDiag(j,5,5)*( lhsK(25,k)
     &                 - BDiag(j,1,5) * lhsK(5,k) 
     &                 - BDiag(j,2,5) * lhsK(10,k) 
     &                 - BDiag(j,3,5) * lhsK(15,k) 
     &                 - BDiag(j,4,5) * lhsK(20,k) )
               enddo
            enddo
c     
      return

      end

      
