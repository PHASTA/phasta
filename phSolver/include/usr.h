/*===========================================================================
 * 
 * "usr.h":  Juin's header funtion
 *
 *===========================================================================
 */

#ifndef	__USR_H__
#define	__USR_H__

#include <FCMangle.h>

/*===========================================================================
 *
 * "UsrHd":  User main struct
 *
 *===========================================================================
 */
typedef struct _Usr {
    int         eqnType ;               /* equation type flag             */
    double*	aperm ;			/* permanent data	          */
    double*	atemp ;			/* temporary data	          */
    double*	resf ;			/* residual vector	          */
    double*	solinc ;		/* solution increment	          */
    double*     flowDiag ;              /* global diagonal of total lhs   */
    double*     sclrDiag ;              /* global diagonal of energy lhs  */
    double*     lesP ;                  /* utility vector for Q = A * P   */
    double*     lesQ ;                  /* utility vector for Q = A * P   */
    int*        iBC  ;                  /* boundary condition code        */
    double*     BC   ;                  /* boundary condition array       */
    int*        iper ;                  /* periodic nodal information     */
    int*        ilwork ;                /* local MPI communication array  */
    int         numpe ;                 /* number of processors           */
    int		nNodes ;		/* number of nodes	          */
    int         nenl ;                  /* number of element nodes        */
    int		nPermDims ;		/* number of permanent data	  */
    int		nTmpDims ;		/* number of temporary data	  */
    int*       	rowp ;		        /* row of p for nonzero's of A    */
    int*       	colm ;		        /* start index for rowp vector    */
    double*     lhsK ;		        /* sparse K matrix (9,nnzeros)    */
    double*     lhsP ;		        /* sparse GoC matrix (4,nnzeros)  */
    double*     lhsS ;
    int*        nnz_tot  ;              /* factor for number of nonzeros) */
    double*     CGsol;                  /* pdot after CG solve */
} Usr ;

typedef struct _Usr* UsrHd ;

/*===========================================================================
 *
 * C to Fortran naming conversion
 *
 *===========================================================================
 */

#define	usrNew		FortranCInterface_GLOBAL_(usrnew,USRNEW)
#define getSol          FortranCInterface_GLOBAL_(getsol,GETSOL)         /* extract Dy from usrHd        */

#define drvlesPrepDiag   FortranCInterface_GLOBAL_(drvlesprepdiag,DRVLESPREPDIAG)
#define drvLesApG        FortranCInterface_GLOBAL_(drvlesapg,DRVLESAPG)
#define drvLesApKG       FortranCInterface_GLOBAL_(drvlesapkg,DRVLESAPKG)
#define drvLesApNGt      FortranCInterface_GLOBAL_(drvlesapngt,DRVLESAPNGT)
#define drvLesApNGtC     FortranCInterface_GLOBAL_(drvlesapngtc,DRVLESAPNGTC)
#define drvLesApFull     FortranCInterface_GLOBAL_(drvlesapfull,DRVLESAPFULL)

#define commOut          FortranCInterface_GLOBAL_(commout,COMMOUT)
#define commIn           FortranCInterface_GLOBAL_(commin,COMMIN)
#define fLesSparseApSclr FortranCInterface_GLOBAL_(flessparseapsclr,FLESSPARSEAPSCLR)
#define fLesSparseApG        FortranCInterface_GLOBAL_(flessparseapg,FLESSPARSEAPG)
#define fLesSparseApKG       FortranCInterface_GLOBAL_(flessparseapkg,FLESSPARSEAPKG)
#define fLesSparseApNGt      FortranCInterface_GLOBAL_(flessparseapngt,FLESSPARSEAPNGT)
#define fLesSparseApNGtC     FortranCInterface_GLOBAL_(flessparseapngtc,FLESSPARSEAPNGTC)
#define fLesSparseApFull     FortranCInterface_GLOBAL_(flessparseapfull,FLESSPARSEAPFULL)

#define drvsclrDiag      FortranCInterface_GLOBAL_(drvsclrdiag,DRVSCLRDIAG)
#define drvLesApSclr     FortranCInterface_GLOBAL_(drvlesapsclr,DRVLESAPSCLR)
#define drvAllreduce     FortranCInterface_GLOBAL_(drvallreduce,DRVALLREDUCE)
#define drvAllreducesclr FortranCInterface_GLOBAL_(drvallreducesclr,DRVALLREDUCESCLR)
#define	flesCp       FortranCInterface_GLOBAL_(flescp,FLESCP)
#define	flesScale    FortranCInterface_GLOBAL_(flesscale,FLESSCALE)
#define	flesScaleCp  FortranCInterface_GLOBAL_(flesscalecp,FLESSCALECP)
#define	flesAdd      FortranCInterface_GLOBAL_(flesadd,FLESADD)
#define	flesSub      FortranCInterface_GLOBAL_(flessub,FLESSUB)
#define	flesDot1     FortranCInterface_GLOBAL_(flesdot1,FLESDOT1)
#define	flesDot2     FortranCInterface_GLOBAL_(flesdot2,FLESDOT2)
#define	flesDaxpy    FortranCInterface_GLOBAL_(flesdaxpy,FLESDAXPY)
#define	flesDxpay    FortranCInterface_GLOBAL_(flesdxpay,FLESDXPAY)
#define	flesInv      FortranCInterface_GLOBAL_(flesinv,FLESINV)
#define	flesZero     FortranCInterface_GLOBAL_(fleszero,FLESZERO)
#define fsclrDiag    FortranCInterface_GLOBAL_(fsclrdiag,FSCLRDIAG)
#define flesApSclr   FortranCInterface_GLOBAL_(flesapsclr,FLESAPSCLR)
#define	fMtxVdimVecMult   FortranCInterface_GLOBAL_(fmtxvdimvecmult,FMTXVDIMVECMULT)
#define	fMtxBlkDot2       FortranCInterface_GLOBAL_(fmtxblkdot2,FMTXBLKDOT2)
#define	fMtxBlkDaxpy      FortranCInterface_GLOBAL_(fmtxblkdaxpy,FMTXBLKDAXPY)
#define	fMtxBlkDyeax      FortranCInterface_GLOBAL_(fmtxblkdyeax,FMTXBLKDYEAX)
#define	fMtxBlkDmaxpy     FortranCInterface_GLOBAL_(fmtxblkdmaxpy,FMTXBLKDMAXPY)
#define	fMtxVdimVecCp     FortranCInterface_GLOBAL_(fmtxvdimveccp,FMTXVDIMVECCP)
#define	fMtxVdimVecDot2   FortranCInterface_GLOBAL_(fmtxvdimvecdot2,FMTXVDIMVECDOT2)
#define	fMtxVdimVecDaxpy  FortranCInterface_GLOBAL_(fmtxvdimvecdaxpy,FMTXVDIMVECDAXPY)
#define ramg_interface    FortranCInterface_GLOBAL_(ramg_interface,RAMG_INTERFACE)
#define ramg_normcheck    FortranCInterface_GLOBAL_(ramg_normcheck,RAMG_NORMCHECK)

/*===========================================================================
 *
 * Function declaration
 *
 *===========================================================================
 */
void 		usrNew(			UsrHd		usrHd,
                                    int*            eqnType,
                                    double*		aperm,
                                    double*		atemp,
                                    double*         resf,
                                    double*         solinc,
                                    double*         flowDiag,
                                    double*         sclrDiag,
                                    double*         lesP,
                                    double*         lesQ,
                                    int*            iBC,
                                    double*         BC,
                                    int*            iper,
                                    int*            ilwork,
                                    int*            numpe,
                                    int*		nNodes,
                                    int*            nenl,
                                    int*		nPermDims,
                                    int*		nTmpDims,
                                    int*		rowp,
                                    int*		colm,
                                    double*		lhsK,
                                    double*		lhsP,
                                    double*         lhsS,
                                    int*             nnz_tot,
                                    double*         CGsol) ;

double*	usrPointer(			UsrHd   usrHd,
                            int	id,
                            int	offset,
                            int	nDims ) ;

void              getSol(      UsrHd            usrHd,
                                        double*          Dy		) ;

/*===========================================================================
 *
 * Fortran Function declaration
 *
 *===========================================================================
 */
double    flesDot1(		double*		a,
                                    int*		m,
                                    int*		n		) ;
double	  flesDot2(		double*		a,
                                    double*		b,
                                    int*		m,
                                    int*		n		) ;

/*===========================================================================
 *
 * End of the file
 *
 *===========================================================================
 */

#endif	/* __USR_H__ */

