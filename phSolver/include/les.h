/****************************************************************************
**  Copyright (c) 1994-2003 ACUSIM Software, Inc.
**  All rights reserved.
**  This source code is confidential and may not be disclosed.
****************************************************************************/

/*===========================================================================
**
** "les.h":  Linear Equation Solvers.
**
** Original: Farzin Shakib (Aug 94)
**===========================================================================
*/

#ifndef	__LES_H__
#define	__LES_H__

/*===========================================================================
 *
 * Get the needed include files
 *
 *===========================================================================
 */
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "usr.h"

/*===========================================================================
 *
 * Main data types
 *
 *===========================================================================
 */
typedef	int	Integer	;		/* integer type			*/
typedef	double	Real ;			/* real    type			*/
typedef	float	Float ;			/* float   type			*/
typedef	char*	String ;		/* string  type			*/
typedef	void*	Data ;			/* data    type			*/
typedef	void	Void ;			/* void    type			*/

#ifdef TRUE
#undef TRUE
#endif
#ifdef FALSE
#undef FALSE
#endif

typedef enum { 
    FALSE = 0 , 
    TRUE  = 1 
} Bool ;				/* boolean type			*/

/*===========================================================================
 *
 * Nil
 *
 *===========================================================================
 */
#define	Nil(T)	((T)NULL)

/*===========================================================================
 *
 * "Les":  Acusim LES data structure
 *
 *===========================================================================
 */
typedef struct _Les* LesHd ;

/*===========================================================================
 *
 * Equation type
 *
 *===========================================================================
 */
#define	LES_EQN_FLOW	1			/* Solve flow vel/pres	*/
#define	LES_EQN_SCLR	2			/* Solve multi-scalar	*/

/*===========================================================================
 *
 * Parameter names
 *
 *===========================================================================
 */
#define	LES_ACT_PRJS		1	/* No. active projection vecs.	*/
#define	LES_ACT_PRES_PRJS	2	/* No. act. pres. prj. vecs.	*/
#define	LES_PRJ_VEC_ID		3	/* Prj vec srcId into permVec	*/
#define	LES_PRES_PRJ_VEC_ID	4	/* Pres. Prj vec srcId 		*/

/*===========================================================================
 *
 * Pointers
 *
 *===========================================================================
 */
#define	LES_RES_PNT		100000		/* residual pointer	*/
#define	LES_SOL_PNT		200000		/* solution pointer	*/

/*===========================================================================
 *
 * Acusim callable funtions in C
 *
 *===========================================================================
 */
LesHd		lesNew(			char*		lmHost,
					Integer		lmPort,
					Integer*	lmNum,
					Integer		eqnType,
					Integer		nDofs,
					Integer		minIters,
					Integer		maxIters,
					Integer		nKvecs,
					Integer		prjFlag,
					Integer		nPrjs,
					Integer		presPrjFlag,
					Integer		nPresPrjs,
					Integer		presPrecFlag,
					Real		tol,
					Real		presTol,
					Integer		verbose,
					Real*		stats,
					Integer*	nPermDims,
					Integer*	nTmpDims	) ;
Void		lesFree(		LesHd		lesHd		) ;
Void		lesSolve(		LesHd		lesHd,
					       UsrHd		usrHd		) ;
Real		lesGetPar(		LesHd		lesHd,
					        Integer		parName		) ;
Void		lesSetPar(		LesHd		lesHd,
					        Integer		parName,
					        Real		parVal		) ;

/*===========================================================================
 *
 * Fortran calling convention
 *
 *===========================================================================
 */
#if	defined(ACUSIM_SGI)     ||	defined(ACUSIM_SGI64)	|| \
	defined(ACUSIM_HAL)	||	defined(ACUSIM_SUN)	|| \
	defined(ACUSIM_ALPHA)	||	defined(ACUSIM_LINUX)	|| \
	defined(ACUSIM_LINUXIPF)
#define	fLesNew		flesnew_
#define	fLesFree	flesfree_
#define	fLesSolve	flessolve_
#elif	defined(ACUSIM_HP) || defined(ibm)
#define	fLesNew		flesnew
#define	fLesFree	flesfree
#define	fLesSolve	flessolve
#elif	defined(ACUSIM_NT)
#define	fLesNew		FLESNEW
#define	fLesFree	FLESFREE
#define	fLesSolve	FLESSOLVE
#else /* dummy */
#define	fLesNew		fLesNewXX
#define	fLesFree	fLesFreeXX
#define	fLesSolve	fLesSolveXX
#endif

/*===========================================================================
 *
 * Acusim callable funtions in Fortran
 *
 *===========================================================================
 */

Void		fLesNew(		Integer*	lesId,
					char*		lmhost,
					Integer*	lmport,
					Integer*	lmNum,
					Integer*	eqnType,
					Integer*	nDofs,
					Integer*	minIters,
					Integer*	maxIters,
					Integer*	nKvecs,
					Integer*	prjFlag,
					Integer*	nPrjs,
					Integer*	presPrjFlag,
					Integer*	nPresPrjs,
					Integer*	presPrecFlag,
					Real*		tol,
					Real*		presTol,
					Integer*	verbose,
					Real*		stats,
					Integer*	nPermDims,
					Integer*	nTmpDims,
					Integer		len_lmHost	) ;
Void		fLesFree(		Integer*	lesId		) ;
Void		fLesSolve(		Integer*	lesId,
					UsrHd		usrHd		) ;

/*===========================================================================
 *
 * Functions to be provided by user
 *
 *===========================================================================
 */
Void		lesPrepDiag(		UsrHd		usrHd		) ;
Void		lesDiagScaleCp(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff,
					Integer		diagOff,
					Integer		nDims		) ;

Void		lesZero(		UsrHd		usrHd,
					Integer		dstId, 
					Integer		nDims		) ;
Void		lesCp(			UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nDims		) ;
Void		lesScale(		UsrHd		usrHd,
					Integer		dstId, 
					Real		coef,
					Integer		nDims		) ;
Void		lesScaleCp(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Real		coef,
					Integer		nDims		) ;
Void		lesAdd(			UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Integer		nDims		) ;
Void		lesSub(			UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Integer		nDims		) ;
Real		lesDot1(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		nDims		) ;
Real		lesDot2(		UsrHd		usrHd,
					Integer		src1Id, 
					Integer		src2Id, 
					Integer		nDims		) ;
Void		lesDaxpy(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Real		coef,
					Integer		nDims		) ;
Void		lesDxpay(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Real		coef,
					Integer		nDims		) ;
Void		lesInv(			UsrHd		usrHd,
					Integer		dstId, 
					Integer		nDims		) ;
Void		lesBlkDot2(		UsrHd		usrHd,
					Integer		src1Id, 
					Integer		src2Id, 
					Real*		values,
					Integer		mDims,
					Integer		nDims		) ;
Void		lesBlkDaxpy(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Real*		coef,
					Integer		mDims,
					Integer		nDims		) ;
Void		lesBlkDyeax(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Real*		coef,
					Integer		mDims,
					Integer		nDims		) ;
Void		lesBlkDmaxpy(		UsrHd		usrHd,
					Integer		srcId, 
					Integer		dstId, 
					Real*		coef,
					Integer		mDims,
					Integer		nDims		) ;
Void		lesVdimCp(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff,
					Integer		nDims		) ;
Void		lesVdimDot2(		UsrHd		usrHd,
					Integer		src1Id,
					Integer		src2Id, 
					Real*		coef,
					Integer		nSrc1Dims,
					Integer		src1Off,
					Integer		nSrc2Dims,
					Integer		src2Off,
					Integer		nDims		) ;
Void		lesVdimDaxpy(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Real*		coef,
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff,
					Integer		nDims		) ;

Void		lesApG(			UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff		) ;
Void		lesApKG(		UsrHd		usrHd,
					Integer		src1Id,
					Integer		src2Id,
					Integer		dstId, 
					Integer		nSrc1Dims,
					Integer		src1Off,
					Integer		nSrc2Dims,
					Integer		src2Off,
					Integer		nDstDims,
					Integer		dstOff		) ;
Void		lesApNGt(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff		) ;
Void		lesApNGtC(		UsrHd		usrHd,
					Integer		src1Id,
					Integer		src2Id,
					Integer		dstId, 
					Integer		nSrc1Dims,
					Integer		src1Off,
					Integer		nSrc2Dims,
					Integer		src2Off,
					Integer		nDstDims,
					Integer		dstOff		) ;
Void		lesApFull(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff		) ;
Void		lesApSclr(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff		) ;
Void	       lesPrecPPE(		UsrHd		usrHd,
					Integer		srcId,
					Integer		dstId, 
					Integer		nSrcDims,
					Integer		srcOff,
					Integer		nDstDims,
					Integer		dstOff		) ;

/*===========================================================================
 *
 * End of the file
 *
 *===========================================================================
 */

#endif	/* __LES_H__ */
