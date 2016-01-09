#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
//MR CHANGE
#include <cstring>
//MR CHANGE END

#include "Input.h"
#include "common_c.h"

#include <FCMangle.h>
#define BC_setVars FortranCInterface_MODULE_(blowercontrol, bc_setvars,		BLOWERCONTROL, BC_SETVARS)

using namespace std; //::cout;
void print_error_code(int ierr);
int SONFATH=0;
extern "C" char phasta_iotype[80];
//extern "C"
extern "C" void BC_setVars(		int*, 		int*, 		
								int*,		int*, 
								double*, 	double*,	double*,
								double*, 	double*,	double*,
								double*, 	double*,	double*, 
								double*		);

int input_fform(phSolver::Input& inp)
{

  int ierr = 0 ;
  int i,j, n_tmp;

  try {
    if(workfc.myrank==workfc.master) {
      printf("\n Complete Filename: %s \n", inp.GetDefaultFileName());
      printf("\n Local Config: %s \n\n", inp.GetUserFileName());
    }

#ifdef AMG
    
	// AMG PARAMETERS 
	
	if ((string)inp.GetValue("Employ AMG") == "True" ) {
	  
	    amgvari.irun_amg = 1;

        amgvari.irun_amg_prec = inp.GetValue("Run AMG As CG-preconditioner");

        amgvarr.strong_eps       = inp.GetValue("Strong Criterion Eps");
        
        amgvarr.ramg_eps         = inp.GetValue("AMG Convergence Eps");
        
        amgvarr.ramg_relax       = inp.GetValue("AMG Relaxation Omega");        
        amgvarr.ramg_trunc       = inp.GetValue("AMG Truncation Set");
        
        amgvari.iamg_verb        = inp.GetValue("AMG Verbosity");
        
        amgvari.iamg_neg_sten    = inp.GetValue("AMG Neg_Sten");           

        amgvari.iamg_nlevel      = inp.GetValue("AMG Nlevel"); 
        
        amgvari.iamg_c_solver  = inp.GetValue("AMG Coarsest Solver");
        
        amgvari.iamg_init = 0;
        amgvari.iamg_setup_frez = inp.GetValue("AMG Freeze Setup");
        if ((string)inp.GetValue("AMG Interpolation Type")=="Standard")
            amgvari.iamg_interp = 1;
        else
            amgvari.iamg_interp = 0;
        amgvari.maxnev = inp.GetValue("AMG GGB nev");
        amgvari.maxncv = inp.GetValue("AMG GGB ncv");
        
        if ((string)inp.GetValue("AMG Smoother Type")=="Gauss-Seidel")
            amgvari.iamg_smoother = 1;
        else if ((string)inp.GetValue("AMG Smoother Type")=="ChebyShev")
            amgvari.iamg_smoother = 2;
        else if ((string)inp.GetValue("AMG Smoother Type")=="MLS")
            amgvari.iamg_smoother = 3;
        amgvarr.ramg_chebyratio       = inp.GetValue("AMG Chebyshev Eigenvalue ratio");

        amgvari.mlsdeg = inp.GetValue("AMG MLS Degree");
        amgvari.iamg_reduce = inp.GetValue("AMG Run Reduced Serial");
	}	
#endif    

/////////////////////////////chen Sep 25 2009  Flow Control Parameters ////////
// Take BC from IC at inlet
      ctrlvari.iI2Binlet = (int)inp.GetValue("Take BC from IC at Inlet");
      if(ctrlvari.iI2Binlet > 0 )
        {
         ctrlvar.inletVelX = inp.GetValue("Inlet Bulk x Velocity");
        }
// set uniform outlet pressure
	//	ctrlvari.isetOutPres = (int)inp.GetValue("Set Outlet Pressure");
	//	if(ctrlvari.isetOutPres > 0 )

    ductvari.isetOutletID  = (int)inp.GetValue("Duct Outlet ID");
    if(ductvari.isetOutletID > 0 )
		ctrlvar.outPres1 = inp.GetValue("Duct Uniform Outlet Pressure");
	
// Override Eddy Vicosity IC and BC
	ductvari.isetEV_IC_BC=(int)inp.GetValue("Override Eddy Viscosity");
	if(ductvari.isetEV_IC_BC==1){
		ductvar.evis_IC_BC=inp.GetValue("Eddy Viscosity Value for Override");
	}

	if(ductvari.isetEVramp = (int)inp.GetValue("Specify Initial Eddy Viscosity Ramp")){
		ductvar.EVrampXmin	=inp.GetValue("Initial Scalar 1 ramp start");
		ductvar.EVrampXmax	=inp.GetValue("Initial Scalar 1 ramp end");
		ductvar.EVrampMax	=inp.GetValue("Initial Scalar 1 high");
		ductvar.EVrampMin   =inp.GetValue("Initial Scalar 1 low");
	}

// set initial condition
	ctrlvari.isetInitial=(int)inp.GetValue("Specify Initial Conditions");
	if(ctrlvari.isetInitial==1){
		ctrlvar.xvel_ini = inp.GetValue("Initial X Velocity");
		ctrlvar.yvel_ini = inp.GetValue("Initial Y Velocity");
		ctrlvar.zvel_ini = inp.GetValue("Initial Z Velocity");
		ctrlvar.temp_ini = inp.GetValue("Initial Temp");
		ctrlvar.pres_ini = inp.GetValue("Initial Pressure");
		ctrlvar.evis_ini = inp.GetValue("Initial Scalar 1");
	}

    

      
//initial condition for duct
	ductvari.isetInitial_Duct=(int)inp.GetValue("Set Initial Condition for Duct");
//inlet condition for duct
	ductvari.isetInlet_Duct=(int)inp.GetValue("Set Inlet Condition for Duct");

    //surfID, t_cycle, t_riseTime, t_fallTime, t_fullOn, vmax, vmin, T, nu, deltaBL, enable
    n_tmp = (int) inp.GetValue("Number of Blower Surfaces");	//BC_setNBlower(&n_tmp);

	if(n_tmp > 0){
		vector<int>    *ivec[3];
		vector<double> *dvec[10];

		for(i = 0; i < 3; i++)	ivec[i] = new vector<int>;
		for(i = 0; i < 10; i++)	dvec[i] = new vector<double>;
		
        *ivec[0] = inp.GetValue("Blower Mode");
		*ivec[1] = inp.GetValue("Blower Surface ID");	
		*ivec[2] = inp.GetValue("Blower Enable");
		
		*dvec[0] = inp.GetValue("Blower Cycle Period");
        *dvec[1] = inp.GetValue("Blower Full On Period");
		*dvec[2] = inp.GetValue("Blower Rise Time");
		*dvec[3] = inp.GetValue("Blower Fall Time");
        *dvec[4] = inp.GetValue("Blower Maximum u_normal");
        *dvec[5] = inp.GetValue("Blower Minimum u_normal");
		*dvec[6] = inp.GetValue("Blower Temperature");	
		*dvec[7] = inp.GetValue("Blower Eddy Viscosity");
		*dvec[8] = inp.GetValue("Blower BL Thickness");		
		*dvec[9] = inp.GetValue("Blower BL Thickness (scalar)");
		
		BC_setVars(	&n_tmp, 		
					&(*ivec[0])[0], 	//mode	
					&(*ivec[1])[0], 	//surfID
					&(*ivec[2])[0],		//enable
					&(*dvec[0])[0],		//t_cycle
					&(*dvec[1])[0],		//t_fullOn
					&(*dvec[2])[0],		//t_riseTime
					&(*dvec[3])[0],		//t_fallTime
				 	&(*dvec[4])[0],		//vmax
					&(*dvec[5])[0],		//vmin
					&(*dvec[6])[0],		//T
					&(*dvec[7])[0],		//nu
				 	&(*dvec[8])[0],		//delta BL velocity
					&(*dvec[9])[0]	);	//delta BL scalar

	}
        
//suction condition for duct
	ductvari.isetSuctionID_Duct=(int)inp.GetValue("Duct Set Suction Surface ID");        //suction patch surface IDs
	if(ductvari.isetSuctionID_Duct > 0){
		ductvari.suctionVbottom     = inp.GetValue("Duct Bottom Suction Normal Velocity");
		ductvari.suctionVside_lower = inp.GetValue("Duct Lower Side Suction Normal Velocity");
		ductvari.suctionVside_upper = inp.GetValue("Duct Upper Side Surface Normal Velocity");
		ductvari.suctionVtop        = inp.GetValue("Duct Top Surface Normal Velocity");	
	}

//  duct geometry type      
	ductvari.iDuctgeometryType = (int)inp.GetValue("Duct Geometry Type");

///////////////////////////////////////////////////////////////////////////////


    // Disabled Features 

    conpar.iALE = inp.GetValue("iALE");
    conpar.icoord = inp.GetValue("icoord");
    conpar.irs = inp.GetValue("irs");
    conpar.iexec = inp.GetValue("iexec");
    timpar.ntseq = inp.GetValue("ntseq");
    solpar.imap = inp.GetValue("imap");


    // Solution Control Keywords

    if((string)inp.GetValue("Equation of State") == "Incompressible") matdat.matflg[0][0] =-1 ;
    if((string)inp.GetValue("Equation of State") == "Compressible") matdat.matflg[0][0] =0;
    inpdat.Delt[0] = inp.GetValue("Time Step Size");
    inpdat.nstep[0] = inp.GetValue("Number of Timesteps");
    if((string)inp.GetValue("Viscous Control")=="Viscous") conpar.navier=1 ; else conpar.navier=0;
   
    if ((string)inp.GetValue("Turbulence Model") == "No-Model" ) {
      turbvari.irans = 0;
      turbvari.iles  = 0;
    } else if ((string)inp.GetValue("Turbulence Model") == "LES" ) {
      turbvari.iles  = 1;
      turbvari.irans = 0;
    } else if ((string)inp.GetValue("Turbulence Model") == "RANS-SA" ) {
      turbvari.iles  = 0;
      turbvari.irans = -1;
    } else if ((string)inp.GetValue("Turbulence Model") == "RANS" ) {
      turbvari.iles  = 0;
      turbvari.irans = -1; // assume S-A for backward compatibility
    } else if ((string)inp.GetValue("Turbulence Model") == "RANS-KE" ) {
      turbvari.iles  = 0;
      turbvari.irans = -2;
    } else if ((string)inp.GetValue("Turbulence Model") == "DES97" ) {
      turbvari.iles  = -1;
      turbvari.irans = -1;
    } else if ((string)inp.GetValue("Turbulence Model") == "DDES" ) {
      turbvari.iles  = -2;
      turbvari.irans = -1;
 
    } else {
      cout << " Turbulence Model: Only Legal Values ( No-Model, LES, RANS-SA, RANS-KE, DES97, DDES )";
      cout << endl;
      exit(1);
    }

 //   if (turbvari.iles*turbvari.irans!=0) turbvar.eles=
 //                                          inp.GetValue("DES Edge Length");

    if (turbvari.irans<0 && turbvari.iles<0)
      turbvar.DES_SA_hmin=(double)inp.GetValue("DES SA Minimum Edge Length");

    int solflow, solheat , solscalr, ilset;
    ((string)inp.GetValue("Solve Flow") == "True")? solflow=1:solflow=0;
    ((string)inp.GetValue("Solve Heat") == "True")? solheat=1:solheat=0;
    //for compressible solheat= False so
    if((string)inp.GetValue("Equation of State") == "Compressible") solheat=0;
    ilset = (int)inp.GetValue("Solve Level Set");
    solscalr = (int)inp.GetValue("Solve Scalars");
    solscalr += ilset;
    if(turbvari.irans == -1) solscalr++;
    if(turbvari.irans == -2) solscalr=solscalr+2;
    if ( solscalr > 4 ) {
      cout << " Only Four Scalars are supported \n";
      cout <<" Please reduce number of scalars \n";
      exit(1);
    }
    inpdat.impl[0] = 10*solflow+solscalr*100+solheat;

    levlset.iLSet = ilset;
    if( ilset > 0) {
    levlset.epsilon_ls = inp.GetValue("Number of Elements Across Interface");
    levlset.epsilon_lsd = inp.GetValue("Number of Elements Across Interface for Redistancing");
    levlset.dtlset = inp.GetValue("Pseudo Time step for Redistancing");
    levlset.iExpLSSclr2 = inp.GetValue("Explicit Solve for Redistance Field");
    levlset.iExpLSSclr1 = inp.GetValue("Explicit Solve for Scalar 1 Field");
    if ((string)inp.GetValue("Apply Volume Constraint") == "True" ) {
      levlset.ivconstraint = 1; } 
    else if((string)inp.GetValue("Apply Volume Constraint") == "False" ) {
      levlset.ivconstraint = 0; }
    else {
      cout << "Apply Volume Constraint: Only Legal Values (True, False) ";
      cout << endl;
      exit(1);
    }   
    }

    vector<double> vec;

    // OUTPUT CONTROL KEY WORDS.

    conpar.necho = inp.GetValue("Verbosity Level");
    outpar.ntout = inp.GetValue("Number of Timesteps between Restarts");
    outpar.nsynciofiles = inp.GetValue("Number of SyncIO Files");
    if((string)inp.GetValue("Print Statistics") == "True") outpar.ioform = 2;
    else outpar.ioform = 1;
  
    if((string)inp.GetValue("Print Wall Fluxes") == "True") outpar.iowflux = 1;
    else outpar.iowflux = 0;

    if((string)inp.GetValue("Print FieldView") == "True") outpar.iofieldv = 1;
    else outpar.iofieldv = 0;

    if((string)inp.GetValue("Print ybar") == "True") outpar.ioybar = 1;
    else outpar.ioybar = 0;

    if((string)inp.GetValue("Print vorticity") == "True") outpar.ivort = 1;
    else outpar.ivort = 0;

    outpar.nstepsincycle = inp.GetValue("Number of Steps in a Cycle");
    outpar.nphasesincycle = inp.GetValue("Number of Phases in a Cycle");
    outpar.ncycles_startphaseavg = inp.GetValue("Number of Initial Cycles to Skip in Phase Average");

    strcpy( outpar.iotype , ((string)inp.GetValue("Data Block Format")).c_str());
    strcpy( phasta_iotype , ((string)inp.GetValue("Data Block Format")).c_str());
    SONFATH = inp.GetValue("Number of Father Nodes");
  
    if((string)inp.GetValue("Print Residual at End of Step") == "True") genpar.lstres = 1;
    else genpar.lstres = 0;
  
    if((string)inp.GetValue("Print Error Indicators") == "True") turbvar.ierrcalc = 1;
    else turbvar.ierrcalc = 0;

    if((string)inp.GetValue("Print Velocity Hessian") == "True") turbvar.ihessian = 1;
    else turbvar.ihessian = 0;

    if ( turbvar.ierrcalc == 1 )
        turbvari.ierrsmooth = inp.GetValue("Number of Error Smoothing Iterations");

    for(i=0;i<MAXSURF+1; i++) aerfrc.nsrflist[i] = 0;
    int nsrfCM = inp.GetValue("Number of Force Surfaces");
    if (nsrfCM > 0) {
      vector<int> ivec = inp.GetValue("Surface ID's for Force Calculation");
      for(i=0; i< nsrfCM; i++){
        aerfrc.nsrflist[ivec[i]] = 1;
        //        cout <<"surface in force list "<< ivec[i] << endl;
      }
      ivec.erase(ivec.begin(),ivec.end());
    }

    aerfrc.isrfIM = inp.GetValue("Surface ID for Integrated Mass");

    outpar.iv_rankpercore = inp.GetValue("Ranks per core");
    outpar.iv_corepernode = inp.GetValue("Cores per node");

    turbvari.iramp=0;
    if((string)inp.GetValue("Ramp Inflow") == "True") turbvari.iramp=1;
    if(turbvari.iramp == 1) {
	vec = inp.GetValue("Mdot Ramp Inflow Start and Stop");
    	for(i=0; i<2 ; i++){
        	turbvar.rampmdot[0][i]=vec[i];
    	}	
    	vec = inp.GetValue("Mdot Ramp Lower FC Start and Stop");
    	for(i=0; i<2 ; i++){
         	turbvar.rampmdot[1][i]=vec[i];
    	}	
    	vec = inp.GetValue("Mdot Ramp Upper FC Start and Stop");
    	for(i=0; i<2 ; i++){
        	turbvar.rampmdot[2][i]=vec[i];
    	}
    }

//Limiting
    vec = inp.GetValue("Limit u1");
    for(i=0; i<3 ; i++){
      turbvar.ylimit[0][i] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

    vec = inp.GetValue("Limit u2");
    for(i=0; i<3 ; i++){
      turbvar.ylimit[1][i] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

    vec = inp.GetValue("Limit u3");
    for(i=0; i<3 ; i++){
      turbvar.ylimit[2][i] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

    vec = inp.GetValue("Limit Pressure");
    for(i=0; i<3 ; i++){
      turbvar.ylimit[3][i] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

    vec = inp.GetValue("Limit Temperature");
    for(i=0; i<3 ; i++){
      turbvar.ylimit[4][i] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

    //Material Properties Keywords 
    matdat.nummat = levlset.iLSet+1;
    if((string)inp.GetValue("Shear Law") == "Constant Viscosity") 
      for(i=0; i < levlset.iLSet+1; i++) matdat.matflg[i][1] = 0;

    if((string)inp.GetValue("Bulk Viscosity Law") == "Constant Bulk Viscosity") 
      for(i=0; i < levlset.iLSet+1; i++) matdat.matflg[i][2] = 0;

    mmatpar.pr = inp.GetValue("Prandtl Number"); 

    if((string)inp.GetValue("Conductivity Law") == "Constant Conductivity") 
      for(i=0; i < levlset.iLSet+1; i++) matdat.matflg[i][3] = 0;

    vec = inp.GetValue("Density");
    for(i=0; i< levlset.iLSet +1 ; i++){
      matdat.datmat[i][0][0] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

    vec = inp.GetValue("Viscosity");
    for(i=0; i< levlset.iLSet +1 ; i++){
      matdat.datmat[i][1][0] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

//      vec = inp.GetValue("Specific Heat");
    for(i=0; i< levlset.iLSet +1 ; i++){
      matdat.datmat[i][2][0] = 0;
    }
//      vec.erase(vec.begin(),vec.end());

    vec = inp.GetValue("Thermal Conductivity");
    for(i=0; i< levlset.iLSet +1 ; i++){
      matdat.datmat[i][3][0] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());
  
    vec = inp.GetValue("Scalar Diffusivity");
    for(i=0; i< solscalr ; i++){
      sclrs.scdiff[i] = vec[i];
    }
    vec.erase(vec.begin(),vec.end());

    if((string)inp.GetValue("Zero Mean Pressure") == "True")
      turbvar.pzero=1;

    turbvar.rmutarget = inp.GetValue("Target Viscosity For Step NSTEP");

    if ( (string)inp.GetValue("Body Force Option") == "None" ) {
      for( i=0; i< levlset.iLSet +1 ; i++)  matdat.matflg[i][4] = 0;
    }
    else if ( (string)inp.GetValue("Body Force Option") == "Vector" ) {
      for( i=0; i< levlset.iLSet +1 ; i++)  matdat.matflg[i][4] = 1;
    }
    else if ( (string)inp.GetValue("Body Force Option") == "User e3source.f" ) {
      for( i=0; i< levlset.iLSet +1 ; i++) matdat.matflg[i][4] = 3;
    }
    else if ( (string)inp.GetValue("Body Force Option") == "Boussinesq" ) {
      for(i=0; i< levlset.iLSet +1 ; i++) matdat.matflg[i][4] = 2;
    }
    else if ( (string)inp.GetValue("Body Force Option") == "Cooling Analytic" ) {
      for(i=0; i< levlset.iLSet +1 ; i++) matdat.matflg[i][4] = 4;
    }
    else if ( (string)inp.GetValue("Body Force Option") == "Cooling Initial Condition" ) {
      for(i=0; i< levlset.iLSet +1 ; i++) matdat.matflg[i][4] = 5;
    }

    // the following block of stuff is common to all cooling type sponges. 
    // Specific stuff belongs in the conditionals above

    if(matdat.matflg[0][4] >=4) {
      spongevar.betamax = inp.GetValue("Maximum Value of Sponge Parameter");
      spongevar.zinsponge = inp.GetValue("Inflow Cooling Sponge Ends at z");
      spongevar.zoutsponge= inp.GetValue("Outflow Cooling Sponge Begins at z");
      spongevar.radsponge = inp.GetValue("Radial Cooling Sponge Begins at r");
      spongevar.grthosponge = inp.GetValue("Sponge Growth Coefficient Outflow");
      spongevar.grthisponge = inp.GetValue("Sponge Growth Coefficient Inflow");


      spongevar.spongecontinuity = 0;
      spongevar.spongemomentum1 = 0;
      spongevar.spongemomentum2 = 0;
      spongevar.spongemomentum3 = 0;
      spongevar.spongeenergy = 0;
 
      if((string)inp.GetValue("Sponge for Continuity Equation") == "True")
	spongevar.spongecontinuity = 1;
      if((string)inp.GetValue("Sponge for x Momentum Equation") == "True")
	spongevar.spongemomentum1 = 1;
      if((string)inp.GetValue("Sponge for y Momentum Equation") == "True")
	spongevar.spongemomentum2 = 1;
      if((string)inp.GetValue("Sponge for z Momentum Equation") == "True")
	spongevar.spongemomentum3 = 1;
      if((string)inp.GetValue("Sponge for Energy Equation") == "True")
	spongevar.spongeenergy = 1;
      
    }

    vec = inp.GetValue("Body Force");
    for(i=0; i< levlset.iLSet +1 ; i++){
      matdat.datmat[i][4][0] = vec[0+i*3];
      matdat.datmat[i][4][1] = vec[1+i*3];
      matdat.datmat[i][4][2] = vec[2+i*3];
    }
    vec.erase(vec.begin(),vec.end());

    vec = inp.GetValue("Body Force Pressure Gradient");
    for(i=0; i< levlset.iLSet +1 ; i++){
      matdat.datmat[i][6][0] = vec[0+i*3];
      matdat.datmat[i][6][1] = vec[1+i*3];
      matdat.datmat[i][6][2] = vec[2+i*3];
    }
    vec.erase(vec.begin(),vec.end());

    if ( (string)inp.GetValue("Surface Tension Option") == "No" ){
        genpar.isurf = 0;
    }
    else if ((string)inp.GetValue("Surface Tension Option") == "Yes" ){
        genpar.isurf = 1;
    }
    else {
      cout << " Surface Tension: Only Legal Values (Yes, No) ";
      cout << endl;
      exit(1);
    }
    if( genpar.isurf > 0) {
      genpar.Bo = inp.GetValue("Bond Number");
    }

    genpar.EntropyPressure = inp.GetValue("Entropy Form of Pressure Constraint on Weight Space");

    
    if ( (string)inp.GetValue("Rotating Frame of Reference") == "True" ) {
      matdat.matflg[0][5] = 1;
      vec = inp.GetValue("Rotating Frame of Reference Rotation Rate");
      matdat.datmat[0][5][0] = vec[0];
      matdat.datmat[0][5][1] = vec[1];
      matdat.datmat[0][5][2] = vec[2];
      vec.erase(vec.begin(),vec.end());
    }
    else {
      matdat.matflg[0][5] = 0;
      matdat.datmat[0][5][0] = 0.;
      matdat.datmat[0][5][1] = 0.;
      matdat.datmat[0][5][2] = 0.;
    }


    //Linear Solver parameters
     conpar.usingpetsc=0;  // default is to have PETSc off
      incomp.iprjFlag = 0; incomp.ipresPrjFlag=0; inpdat.svLSFlag=0;
      inpdat.leslib=0;
    if( (string)inp.GetValue("Solver Type") =="svLS" ){
      inpdat.svLSFlag = 1; }
    if( (string)inp.GetValue("Solver Type") =="ACUSIM with P Projection" ){
      inpdat.leslib=1; incomp.iprjFlag = 0; incomp.ipresPrjFlag=1;}
    else if ( (string)inp.GetValue("Solver Type") =="ACUSIM" ){
      inpdat.leslib=1; incomp.iprjFlag = 0; incomp.ipresPrjFlag=0;}
    else if( (string)inp.GetValue("Solver Type") =="ACUSIM with Velocity Projection" ){
      inpdat.leslib=1; incomp.iprjFlag = 1; incomp.ipresPrjFlag=0;}
    else if( (string)inp.GetValue("Solver Type") =="ACUSIM with Full Projection" ){
      inpdat.leslib=1; incomp.iprjFlag = 1; incomp.ipresPrjFlag=1;}
    else if( (string)inp.GetValue("Solver Type") =="GMRES Matrix Free"){ 
      inpdat.impl[0] += 10*solflow;}
    else if( (string)inp.GetValue("Solver Type") =="GMRES EBE"){ 
      inpdat.impl[0] += 20*solflow;}
    else if( (string)inp.GetValue("Solver Type") =="PETSc"){ 
      conpar.usingpetsc=1;}
    //GMRES sparse is assumed default and has the value of 10, MFG 20,
    // EBE 30


    //    inpdat.niter[0] = inp.GetValue("Number of Solves per Time Step");
    solpar.nGMRES = inp.GetValue("Number of GMRES Sweeps per Solve");
    solpar.Kspace = inp.GetValue("Number of Krylov Vectors per GMRES Sweep");
    inpdat.LHSupd[0] = inp.GetValue("Number of Solves per Left-hand-side Formation");
    inpdat.epstol[0] = inp.GetValue("Tolerance on Momentum Equations");
    incomp.prestol = inp.GetValue("Tolerance on ACUSIM Pressure Projection"); 
    incomp.minIters = inp.GetValue("Minimum Number of Iterations per Nonlinear Iteration");
    incomp.maxIters = inp.GetValue("Maximum Number of Iterations per Nonlinear Iteration");
    inpdat.deltol[0][0]=inp.GetValue("Velocity Delta Ratio"); 
    inpdat.deltol[1][0]=inp.GetValue("Pressure Delta Ratio"); 
    incomp.nPrjs = inp.GetValue("Number of Velocity Projection Vectors");
    incomp.nPresPrjs = inp.GetValue("Number of Pressure Projection Vectors");
    incomp.iverbose = inp.GetValue("ACUSIM Verbosity Level");

    if(solheat==1){ 
      inpdat.epstol[1]=inp.GetValue("Temperature Solver Tolerance");
      inpdat.LHSupd[1]=inp.GetValue("Number of Solves of Temperature per Left-hand-side Formation");
    }

    // The following is where you should put any inputs that are able to 
    // input differently for each scalar.  It is a little tedious in the code 
    // but it should make the solver.inp easier to understand. Note this will 
    // require some care with regression tests.


    if(solscalr>0){
      inpdat.epstol[2]=inp.GetValue("Scalar 1 Solver Tolerance");
      inpdat.LHSupd[2]=inp.GetValue("Number of Solves of Scalar 1 per Left-hand-side Formation");

      vec = inp.GetValue("Limit Scalar 1");
      for(i=0; i<3 ; i++){
        turbvar.ylimit[5][i] = vec[i];
      }
      vec.erase(vec.begin(),vec.end());
    } 

    if(solscalr>1){
      inpdat.epstol[3]=inp.GetValue("Scalar 2 Solver Tolerance");
      inpdat.LHSupd[3]=inp.GetValue("Number of Solves of Scalar 2 per Left-hand-side Formation");

      vec = inp.GetValue("Limit Scalar 2");
      for(i=0; i<3 ; i++){
        turbvar.ylimit[6][i] = vec[i];
      }
      vec.erase(vec.begin(),vec.end());
    } 

    if(solscalr>2){
      inpdat.epstol[4]=inp.GetValue("Scalar 3 Solver Tolerance");
      inpdat.LHSupd[4]=inp.GetValue("Number of Solves of Scalar 3 per Left-hand-side Formation");

      vec = inp.GetValue("Limit Scalar 3");
      for(i=0; i<3 ; i++){
        turbvar.ylimit[7][i] = vec[i];
      }
      vec.erase(vec.begin(),vec.end());
    } 

    if(solscalr>3){
      inpdat.epstol[5]=inp.GetValue("Scalar 4 Solver Tolerance");
      inpdat.LHSupd[5]=inp.GetValue("Number of Solves of Scalar 4 per Left-hand-side Formation");

      vec = inp.GetValue("Limit Scalar 4");
      for(i=0; i<3 ; i++){
        turbvar.ylimit[8][i] = vec[i];
      }
      vec.erase(vec.begin(),vec.end());
    }
 
    // DISCRETIZATION CONTROL
    
    genpar.ipord = inp.GetValue("Basis Function Order");
    if((string)inp.GetValue("Time Integration Rule") == "First Order")
      inpdat.rhoinf[0] = -1 ;
    else inpdat.rhoinf[0] = (double)inp.GetValue("Time Integration Rho Infinity");
    if((string)inp.GetValue("Predictor at Start of Step")=="Same Velocity")
      genpar.ipred = 1;
    if((string)inp.GetValue("Predictor at Start of Step")=="Zero Acceleration")
      genpar.ipred = 2;
    if((string)inp.GetValue("Predictor at Start of Step")=="Same Acceleration")
      genpar.ipred = 3;
    if((string)inp.GetValue("Predictor at Start of Step")=="Same Delta")
      genpar.ipred = 4;
    
    if((string)inp.GetValue("Weak Form") == "Galerkin")
      solpar.ivart = 1;
    if((string)inp.GetValue("Weak Form") == "SUPG")
      solpar.ivart = 2;

    if((string)inp.GetValue("Flow Advection Form") == "Convective")
      solpar.iconvflow = 2;
    else if((string)inp.GetValue("Flow Advection Form") == "Conservative")
      solpar.iconvflow = 1;
    if((string)inp.GetValue("Scalar Advection Form") == "Convective")
      solpar.iconvsclr = 2;
    else if((string)inp.GetValue("Scalar Advection Form") == "Conservative")
      solpar.iconvsclr = 1;
    if((string)inp.GetValue("Use Conservative Scalar Convection Velocity") == "True")
      sclrs.consrv_sclr_conv_vel = 1;
    else if((string)inp.GetValue("Use Conservative Scalar Convection Velocity") == "False")
      sclrs.consrv_sclr_conv_vel = 0;
    // TAU INPUT 
    if((string)inp.GetValue("Tau Matrix") == "Diagonal-Shakib")
      genpar.itau = 0;
    else  if((string)inp.GetValue("Tau Matrix") == "Diagonal-Franca")
      genpar.itau =1;
    else if((string)inp.GetValue("Tau Matrix") == "Diagonal-Jansen(dev)") 
      genpar.itau = 2;
    else if((string)inp.GetValue("Tau Matrix") == "Diagonal-Compressible")
      genpar.itau = 3;
    else if((string)inp.GetValue("Tau Matrix") == "Matrix-Mallet") 
      genpar.itau = 10;
    else if((string)inp.GetValue("Tau Matrix") == "Matrix-Modal")
      genpar.itau = 11;

    genpar.dtsfct = inp.GetValue("Tau Time Constant");
    genpar.taucfct = inp.GetValue("Tau C Scale Factor");

	genpar.iLHScond = inp.GetValue("LHS BC heat flux enable");

    // FLOW DISCONTINUITY CAPTURING

      if((string)inp.GetValue("Discontinuity Capturing") == "Off") solpar.iDC = 0;
    else if((string)inp.GetValue("Discontinuity Capturing") == "DC-mallet") solpar.iDC = 1;
    else if((string)inp.GetValue("Discontinuity Capturing") == "DC-quadratic") solpar.iDC = 2;
   else if((string)inp.GetValue("Discontinuity Capturing") == "DC-minimum") solpar.iDC = 3;    
    else {
      cout<< "Condition not defined for Discontinuity Capturing \n ";
      exit(1);
    }

    // SCALAR DISCONTINUITY CAPTURING

      vector<int> ivec = inp.GetValue("Scalar Discontinuity Capturing");
      for(i=0; i< 2; i++)  solpar.idcsclr[i] = ivec[i];
      ivec.erase(ivec.begin(),ivec.end());
 

//        if((string)inp.GetValue("Scalar Discontinuity Capturing") == "No") solpar.idcsclr = 0;
//      else if((string)inp.GetValue("Scalar Discontinuity Capturing") == "1") solpar.idcsclr = 1; 
//   else if((string)inp.GetValue("Scalar Discontinuity Capturing") == "2") solpar.idcsclr = 2; 
//   else {
//        cout<< "Condition not defined for Scalar Discontinuity Capturing \n ";
//        exit(1);
//      }     
    if((string)inp.GetValue("Include Viscous Correction in Stabilization") == "True")
      {  
        if(genpar.ipord == 1 ) genpar.idiff = 1;
        else genpar.idiff = 2;
      }
    else { genpar.idiff = 0;}

// ------ Only For duct S duct Project ------------------------------------------------
    genpar.irampViscOutlet = (int)inp.GetValue("Ramp Up Viscosity Near Outlet");

    genpar.istretchOutlet = (int)inp.GetValue("Stretch X Coordinate Near Outlet");
// -----------------------------------------------------------------------------------

    genpar.iremoveStabTimeTerm = (int)inp.GetValue("Remove Time Term from Stabilization");

    timdat.flmpl = inp.GetValue("Lumped Mass Fraction on Left-hand-side");
    timdat.flmpr = inp.GetValue("Lumped Mass Fraction on Right-hand-side");

    timdat.iCFLworst = 0;
    if((string)inp.GetValue("Dump CFL") == "True")
      timdat.iCFLworst = 1;

    intdat.intg[0][0]=inp.GetValue("Quadrature Rule on Interior");
    intdat.intg[0][1]=inp.GetValue("Quadrature Rule on Boundary");
    genpar.ibksiz = inp.GetValue("Number of Elements Per Block");

    ((string)inp.GetValue("Turn Off Source Terms for Scalars") 
         == "True") ? sclrs.nosource = 1 : sclrs.nosource = 0;
    sclrs.tdecay=inp.GetValue("Decay Multiplier for Scalars");

    // TURBULENCE MODELING PARAMETER
    int tpturb = turbvari.iles-turbvari.irans;
    int ifrule;
    if( tpturb != 0 ){


      turbvari.nohomog =inp.GetValue("Number of Homogenous Directions");

      if((string)inp.GetValue("Turbulence Wall Model Type") == "Slip Velocity") turbvar.itwmod = 1;
      else if((string)inp.GetValue("Turbulence Wall Model Type") == "Effective Viscosity") turbvar.itwmod = 2; 
      else  turbvar.itwmod = 0;
      if (turbvari.irans < 0) turbvar.itwmod = turbvar.itwmod*(-1);
      ifrule  = inp.GetValue("Velocity Averaging Steps");
      turbvar.wtavei =(ifrule >0)? 1.0/ifrule : -1.0/ifrule;
 
      if(turbvari.iles == 1){
        
        if((string)inp.GetValue("Dynamic Model Type") == "Bardina") turbvari.iles += 10;
        else if((string)inp.GetValue("Dynamic Model Type") == "Projection") turbvari.iles += 20;

        ifrule = inp.GetValue("Filter Integration Rule");
        turbvari.iles += ifrule-1;
        ifrule = inp.GetValue("Dynamic Model Averaging Steps");
        turbvar.dtavei = (ifrule >0)? 1.0/ifrule : -1.0/ifrule;
        turbvar.fwr1 = inp.GetValue("Filter Width Ratio");
        turbvar.flump = inp.GetValue("Lumping Factor for Filter");


        if ((string)inp.GetValue("Model Statistics") == "True" ) {
          turbvari.modlstats = 1; } 
        else {
          turbvari.modlstats = 0; }   
 
        if ((string)inp.GetValue("Double Filter") == "True" ) {
          turbvari.i2filt = 1; } 
        else {
          turbvari.i2filt = 0; }  

        if ((string)inp.GetValue("Model/SUPG Dissipation") == "True" ) {
          turbvari.idis = 1; } 
        else {
          turbvari.idis = 0; }


        if((string)inp.GetValue("Dynamic Model Type") == "Standard") {

          if((string)inp.GetValue("Dynamic Sub-Model Type") == "None") 
            turbvari.isubmod = 0;
          else if((string)inp.GetValue("Dynamic Sub-Model Type") =="DFWR")
            turbvari.isubmod = 1;
          else if((string)inp.GetValue("Dynamic Sub-Model Type") =="SUPG") 
            turbvari.isubmod = 2;
        }
        else if((string)inp.GetValue("Dynamic Model Type") == "Projection") {

          if((string)inp.GetValue("Projection Filter Type") == "Linear")
            turbvari.ifproj = 0;
          else if((string)inp.GetValue("Projection Filter Type") =="Quadratic") 
            turbvari.ifproj = 1;          

          if((string)inp.GetValue("Dynamic Sub-Model Type") == "None")
            turbvari.isubmod = 0;
          else if((string)inp.GetValue("Dynamic Sub-Model Type") =="ConsistentProj") 
            turbvari.isubmod = 1;          
        }

      }
    }
  
    // SPEBC MODELING PARAMETERS

    if ( (spebcvr.irscale = inp.GetValue("SPEBC Model Active")) >= 0 ){

      ifrule  = inp.GetValue("Velocity Averaging Steps");
      turbvar.wtavei =(ifrule >0)? 1.0/ifrule : 1.0/inpdat.nstep[0];
      spebcvr.intpres = inp.GetValue("Interpolate Pressure");
      spebcvr.plandist = inp.GetValue("Distance between Planes");
      spebcvr.thetag  = inp.GetValue("Theta Angle of Arc");
      spebcvr.ds = inp.GetValue("Distance for Velocity Averaging");
      spebcvr.tolerence = inp.GetValue("SPEBC Cylindrical Tolerance");
      spebcvr.radcyl = inp.GetValue("Radius of recycle plane");
      spebcvr.rbltin  = inp.GetValue("Inlet Boundary Layer Thickness");
      spebcvr.rvscal  = inp.GetValue("Vertical Velocity Scale Factor");
    } 
      
    // CARDIOVASCULAR MODELING PARAMETERS
    if ( (string)inp.GetValue("Time Varying Boundary Conditions From File") == "True") 
      nomodule.itvn = 1;
    else 
      nomodule.itvn = 0;

    if ( nomodule.itvn ==1)
      nomodule.bcttimescale = inp.GetValue("BCT Time Scale Factor");

    nomodule.ipvsq=0;
    if( (nomodule.icardio = inp.GetValue("Number of Coupled Surfaces")) ){
      if ( nomodule.icardio > MAXSURF ) {
        cout << "Number of Coupled Surfaces > MAXSURF \n";
        exit(1);
      } 
      if ( (string)inp.GetValue("Pressure Coupling") == "None") 
        nomodule.ipvsq=0;
      if ( (string)inp.GetValue("Pressure Coupling") == "Explicit") 
        nomodule.ipvsq=1;
      if ( (string)inp.GetValue("Pressure Coupling") == "Implicit") 
        nomodule.ipvsq=2;
      if ( (string)inp.GetValue("Pressure Coupling") == "P-Implicit") 
        nomodule.ipvsq=3;

      if( (nomodule.numResistSrfs=inp.GetValue("Number of Resistance Surfaces")) ){
          ivec = inp.GetValue("List of Resistance Surfaces");          
          for(i=0;i<MAXSURF+1; i++) nomodule.nsrflistResist[i] = 0;
          for(i=0; i< nomodule.numResistSrfs; i++){
              nomodule.nsrflistResist[i+1] = ivec[i];
          }
          vec = inp.GetValue("Resistance Values");
          for(i =0; i< MAXSURF+1 ; i++) nomodule.ValueListResist[i] = 0;
          for(i =0; i< nomodule.numResistSrfs ; i++) nomodule.ValueListResist[i+1] = vec[i];
          vec.erase(vec.begin(),vec.end());
      }
      if( (nomodule.numImpSrfs=inp.GetValue("Number of Impedance Surfaces")) ){
          ivec = inp.GetValue("List of Impedance Surfaces");
          for(i=0;i<MAXSURF+1; i++) nomodule.nsrflistImp[i] = 0;
          for(i=0; i< nomodule.numImpSrfs; i++){
              nomodule.nsrflistImp[i+1] = ivec[i];
          }
          if ( (string)inp.GetValue("Impedance From File") == "True")
              nomodule.impfile = 1; else nomodule.impfile = 0;
      }
      if( (nomodule.numRCRSrfs=inp.GetValue("Number of RCR Surfaces")) ){
          ivec = inp.GetValue("List of RCR Surfaces");
          for(i=0;i<MAXSURF+1; i++) nomodule.nsrflistRCR[i] = 0;
          for(i=0; i< nomodule.numRCRSrfs; i++){
              nomodule.nsrflistRCR[i+1] = ivec[i];
          }
          if ( (string)inp.GetValue("RCR Values From File") == "True")
              nomodule.ircrfile = 1; else nomodule.ircrfile = 0;
      }
    }
    nomodule.ideformwall = 0;
    if((string)inp.GetValue("Deformable Wall")=="True"){
        nomodule.ideformwall = 1;
        nomodule.rhovw = inp.GetValue("Density of Vessel Wall");
        nomodule.thicknessvw = inp.GetValue("Thickness of Vessel Wall");
        nomodule.evw = inp.GetValue("Young Mod of Vessel Wall");
        nomodule.rnuvw = inp.GetValue("Poisson Ratio of Vessel Wall");
        nomodule.rshearconstantvw = inp.GetValue("Shear Constant of Vessel Wall");
        if((string)inp.GetValue("Wall Mass Matrix for LHS") == "True") nomodule.iwallmassfactor = 1;
        else nomodule.iwallmassfactor = 0;
        if((string)inp.GetValue("Wall Stiffness Matrix for LHS") == "True") nomodule.iwallstiffactor = 1;
        else nomodule.iwallstiffactor = 0; 
    }
    nomodule.iviscflux = 1;
    if((string)inp.GetValue("Viscous Flux Flag") == "True") nomodule.iviscflux = 1;
    if((string)inp.GetValue("Viscous Flux Flag") == "False") nomodule.iviscflux = 0;

   
    // Scaling Parameters Keywords

    outpar.ro = inp.GetValue("Density");
    outpar.vel = inp.GetValue("Velocity");
    outpar.press = inp.GetValue("Pressure");
    outpar.temper = inp.GetValue("Temperature");
    outpar.entrop = inp.GetValue("Entropy");

    // Step Sequencing
 

    ivec = inp.GetValue("Step Construction");
    sequence.seqsize = ivec.size();
    if( sequence.seqsize > 100 || sequence.seqsize < 2 )
     cerr<<"Sequence size must be between 2 and 100 "<<endl;
   
    for(i=0; i< sequence.seqsize; i++)
      sequence.stepseq[i] = ivec[i];

  }
  catch ( exception &e ) {
    cout << endl << "Input exception: " << e.what() << endl << endl;
    ierr = 001;
    print_error_code(ierr);
    return ierr;
  }

  return ierr;
  
}

void print_error_code(int ierr) {
  /*
    Return Error codes:
    0xx         Input error
    1xx         Solution Control
    105         Turbulence Model not supported

    2xx         Material Properties

    3xx         Output Control

    4xx         Discretization Control

    5xx         Scaling Parameters

    6xx         Linear Solver Control
    601         linear solver type not supported
  */
  cout << endl << endl << "Input error detected: " << endl << endl;
  if ( ierr == 001 ) {
    cout << endl << "Input Directive not understood" << endl << endl;
  }
  if ( ierr == 105 ) {
    cout << endl << "Turbulence Model Not Supported" << endl << endl;
  }
  if ( ierr == 601 ) {
    cout << endl << "Linear Solver Type Not Supported" << endl << endl;
  }

}
