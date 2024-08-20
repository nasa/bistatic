/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                     //
// Property of National Aeronautics and Space Administration.                                          //
//                                                                                                     //
// National Aeronautics and Space Administration CONFIDENTIAL                                          //
//                                                                                                     // 
// NOTICE:  All information contained herein is, and remains                                           //
// the property of National Aeronautics and Space Administration SAC and its approved contractors. The //
// intellectual and technical concepts contained herein are proprietary to National Aeronautics and    //
// Space Administration.  Dissemination of this information or reproduction of this material           //
// is strictly forbidden unless prior written permission is obtained from National Aeronautics and     // 
// Space Administration.                                                                               //
//                                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                     //
// Function Inputs: A [Nx3] and corresponding vectors in matrix B [Nx3]                                //
//                                                                                                     //
// Function Outputs: THETA = [Nx1] (deg)                                                               //
//                                                                                                     //
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                     //
//   Function Description                                                                              //
//   Compute angles between vectors in matrix A [Nx3] and corresponding vectors in matrix B [Nx3] in   //
//   which in both matrices the vector components are arranged row-wise [X Y Z]. Return angles between //
//   corresponding vectors in vector set A and vector set B in the [Nx1] vector THETA (deg)            //
//   Note that A and B MUST BE defined in the SAME coordinate frame (e.g. ECEF, ECI, etc) and assumed  //
//   to have the same origin                                                                           //
//       																							   //
//   parallelized version of VECANG.m                                                                  //
//                                                                                                     //
// Last Edit: $Date$                                                                                   //
// ID: $Id$                                                                                            //
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// THETA = VECANGmexmult(A,B)

#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <strsafe.h>
#include <math.h>
#include <emmintrin.h>
#include <errno.h>

#include "mex.h"
#define MAX_THREADS 8
#define block 1000 // each thread handles this many rows

HANDLE ghMutex; 
DWORD ThreadProc (LPVOID lpdwThreadParam );

//% This code calculate bi-static scattering from rough surface
//% using the two scale scattering model
//% Input Parameters are:
//%  
//%  epsr  : complex relative permittivity
//%  U     : Wind speed at 10m (m/s)
//%  k     :  RF wavenumber (1/m)
//%  Q     :  Inverse wave age
//%  mu2    : the upwind mean square slope.
//%  mc2    : the crosswind mean square slopes. 
//%  thetai : zenith incident angle (Deg)
//%  phi    : azimuth incident angle (Deg)
//%  thetas : scattering angle (Deg)
//%  phs    : azimuth incident angle (Deg)
//%  psi    : The orientation of upwing direction with respect to x axis
//%           which is taken to equal to phi (Deg)
//%
//[Gvv,Gvh,Ghv,Ghh]=sea_sur_short(epsr,U,Q,k,mu2,mc2,kd,thi,ths,phi,phs,psi)
double *Gvv, *Gvh, *Ghv, *Ghh,*UPtr,*QPtr,*kPtr,*mu2Ptr,*mc2Ptr,*kdPtr,*thi,*ths,*phi,*phs,*psi;
mxArray * epsr;

double *xu,*wu,*xc,*wc,scmax,sumax,scmin,DK,SK;
int countThread[MAX_THREADS+1];
size_t numRows;

long CurrentIndex;
//Global variable Shared by all threads
int nGlobalCount = 0;
//Main function which starts out each thread

void complex_sqrt(double *a, double *b, double c, double d){
    *a = sqrt((c+sqrt(c*c+d*d))/2.0);
    *b = d/fabs(d) * sqrt((-c+sqrt(c*c+d*d))/2.0); 
}

double complex_abs(double c, double d){
    return sqrt(c*c+d*d);
}

//[ghh,gvh,ghv,gvv]=small_pert_model(epsr,thi,ths,phi,phs)
void small_pert_model(double *ghhR,double *ghhI,double *gvhR,double *gvhI,double *ghvR,double *ghvI,double *gvvR,double *gvvI, mxArray *epsrVal, double thi,double ths,double phi,double phs){
    
    double deg = 57.295779513082320;
    double Rad = 0.017453292519943;
    double thi_val, ths_val, phi_val, phs_val;
    double ci, si, ss,cs,phsi,sphi,cphi,Di,D_R,D_I,DsR,DsI,DiR,DiI,xR,xI,yR,yI;
    double *EsprR, *EsprI,EsprR_Minus_1,ghhIntermediateR,ghhIntermediateI;

    thi_val = thi;
    ths_val = ths;
    phi_val = phi;
    phs_val = phs;   
    //mxComplexDouble * epsrComplex = mxGetComplexDoubles(epsrVal);
    EsprR = mxGetPr(epsrVal);
    EsprI = mxGetPi(epsrVal);

    /*
    thi_val=thi_val*Rad;
    ths_val=ths_val*Rad;
    phi_val=phi_val*Rad;
    phs_val=phs_val*Rad;
    */
    
    ci=cos(thi_val);
    si=sin(thi_val);
    ss=sin(ths_val);
    cs=cos(ths_val);

    phsi = phs_val-phi_val;
    sphi=sin(phsi);
    cphi=cos(phsi);
    
    // *ghh, *gvh, *ghv, *gvv
    complex_sqrt(&DsR,&DsI,*EsprR-ss*ss,*EsprI);
    complex_sqrt(&DiR,&DiI,*EsprR-si*si,*EsprI);
    
    //D=(cs+Ds)*(ci+Di); 
    D_R=(cs+DsR)*(ci+DiR)-(DsI*DiI);
    D_I=(cs+DsR)*DiI + (ci+DiR)*DsI;
    
    EsprR_Minus_1 = *EsprR-1;
    
    //*(ghh)=(epsr_val-1)*cphi/D;   // hh
    ghhIntermediateR = EsprR_Minus_1*cphi;// multiple by real cphi
    ghhIntermediateI = *EsprI * cphi;
    //compute real  (re *re - imag*conj_imag)/(D_re*D_re+D_imag*D_imag)
    *ghhR =  (ghhIntermediateR * D_R - ghhIntermediateI * (-D_I))/(D_R * D_R + D_I *D_I);
    *ghhI = (ghhIntermediateR * (-D_I) + ghhIntermediateI * (D_R))/(D_R * D_R + D_I *D_I);
    //*ghhR=D_R;
    //*ghhI=D_I;
    
    //D = (epsr*cs + Ds)*(ci+Di); 
    D_R = (*EsprR * cs + DsR)*(ci+DiR)-(*EsprI*cs+DsI)*(DiI);
    D_I = (*EsprI * cs + DsI)*(ci+DiR)+(*EsprR*cs+DsR)*(DiI);  
    
    //*(gvh)=-(epsr_val-1)*Ds*sphi/D;  // vh
    xR = -sphi*(EsprR_Minus_1 * DsR-*EsprI * DsI);
    xI = -sphi*(EsprR_Minus_1 * DsI + *EsprI * DsR);
    *(gvhR) = (xR*D_R + xI*(D_I))/(D_R*D_R+D_I*D_I);  //it is + xI*D_I to get conj
    *(gvhI) = (-xR*D_I + xI*(D_R))/(D_R*D_R+D_I*D_I);  //it is + xI*D_I to get conj

    //D = (cs + Ds)*(epsr*ci+Di);
    D_R = ((cs+DsR)* (*EsprR*ci+DiR))-((DsI) * ( *EsprI * ci + DiI));
    D_I = (DsI)*((*EsprR)*ci+DiR)+((cs+DsR)*(*EsprI * ci + DiI));
 
    //ghv=(epsr-1)*Di*sphi/D;   %hv
     xR = sphi*(EsprR_Minus_1*DiR-*EsprI*DiI); //complex multiply DiR and Espr-1
     xI = sphi*(EsprR_Minus_1*DiI+ *EsprI*DiR);   
    *(ghvR) = (xR*D_R + xI*D_I)/(D_R*D_R+D_I*D_I);  //it is + xI*D_I to get conj
    *(ghvI) = (-xR*D_I + xI*D_R)/(D_R*D_R+D_I*D_I);  //it is + xI*D_I to get conj 
   
    //D = (epsr*cs + Ds)*(epsr*ci + Di);
    D_R = (*EsprR*cs +DsR)*(*EsprR*ci+DiR) - (*EsprI*cs +DsI)*(*EsprI*ci+DiI) ;
    D_I = (*EsprR*cs +DsR)*(*EsprI*ci+DiI) + (*EsprI*cs+DsI)*(*EsprR*ci+DiR) ;
    
    //*(gvv)=(epsr_val-1)*(epsr*si*ss-Ds*Di*cphi)/D; // vv   
    //compute x = (epsr*si*ss-Ds*Di*cphi)
     xR = *EsprR * si *ss - cphi * (DsR*DiR - DsI*DiI);
     xI = *EsprI * si *ss - cphi * (DsR*DiI + DsI*DiR);
    //compute y = x/D
    yR = (xR*D_R+xI*D_I)/(D_R*D_R +D_I*D_I);
    yI = (xR*-D_I+xI*D_R)/(D_R*D_R +D_I*D_I);
     
     //compute epsr_val-1 * Y
    *(gvvR) = EsprR_Minus_1*yR - *EsprI * yI;  
    *(gvvI) = EsprR_Minus_1*yI + *EsprI * yR;
}

//% sea surface spectra 
void sea_sur_spectra(double *S,double*DK,double U,double k,double Q){
    // This code calculate sea surface spectra as a function of wavenumber k
    // for several wind speed values at 10 m
    // Based on Eqs. (D-3)  - (D -9)
    //format long
    // Input:
    //  U   wind speed at 10 m (m/s)
    //  k   spatial wavenumber
    //  Q   inverse wave age
    //
    // Output
    // S    value of isotropic (omni-directional) spectrum function
    // DK   amplitude of the directional angular part of the spreading function
    double g=9.81;
    double km=364.52;
    double xq,u,kp,am,CK,Bx,Bh,g1,g2,G,Snm,yk,gama,k3;
    xq=Q/sqrt(10);
    
    u=U*sqrt(0.001*(0.81+0.065*U));
    kp=g*(Q/U)*(Q/U);
    am=0.014*u/0.232;
    //
    CK=sqrt(g/k*(1 + pow((k/km),2.0)));
    //
    Bx=0.003*sqrt(Q)*U/(Q*CK);
    Bx=Bx*exp(-xq*(sqrt(k/kp)-1));
    //
    Bh=0.5*am*(0.232/CK);
    Bh=Bh*exp(-0.25*pow((k/km-1),2.0));
     //
    g1=pow((sqrt(k/kp)-1),2.0);
    if(Q<5){
      g2=2*pow(0.08*(1+4/pow(Q,3.0)),2.0);
    }
    else{
      g2=0.0512;
    }
    gama=exp(-g1/g2);
    
    if(Q<1){// %RM fix
        G=1.7;
    }
    else if(Q<5){
        G=1.7 + 6.0*log(Q);
    }
    else{
        G=2.7*pow(Q,0.57);
    }
    Snm= (Bx+Bh)*exp(-1.25*(kp/k)*(kp/k))*pow(G,gama);
    k3=k*k*k;
    *S=Snm/k3;
    ///
    yk=log(2)/4 + 4*pow((Q*CK/U),2.5)+(0.13*u/0.232)*pow((0.232/CK),2.5);
    *DK=tanh(yk);
}

//function [x w]=gauss(n)
void gauss(int n,double *x,double *w){
    double eps = 0.1E-9;
    double nmax = 10;
    double ip   = 512;
    double *c1,*p,*c2;
    double xnow,pdir;
    double pi4 = 0.78539816339744830961566084581988;
    double pi = 3.1415926535897932384626433832795;
    double constant;
    int k,n1,n2,num_root,ncount;
    double xinc;
    int i;
    bool found,a,b;
    n1 = n+1;
    c1 = mxCalloc(n, sizeof(double));
    c2 = mxCalloc(n, sizeof(double));
    p = mxCalloc(n1, sizeof(double));
    
    *p = 1;

    // Calculate the coefficients for the Legendre Poly.recursive formula
    *(c1) = 0;
    *(c2) = 0;
    for(int ii=2; ii<=n;ii++){
        *(c1+ii-1)=(double)(2*ii-1)/(double)ii;
        *(c2+ii-1)=(double)(ii-1)/(double)ii;
    }
    
   //mexPrintf("return\n");
 
    // Initial Constants
 
    // pi4 = pi/4;
    constant = 1/((double)n+0.5);

    // Determine number of roots (num_root) needed to be calculated
     n2 = n>>1;
     if(2*n2 == n){
         num_root = n2;
     }
     else{
         num_root = n2 + 1;
     }
    
    // Main loop begins here
    // mexPrintf("num_root %d\n", num_root);
     for(i=0; i<num_root;i++){
         k = n-i-1;
         ncount = 0;
         xinc = 1;
    
    // Use Newton's method and a good initial guess to find root

         xnow = cos((((double)i+1.0)*pi-pi4)*constant);
         found = 0;
         
         while (!found){   
             ncount++;
             p[1] = xnow;
             //mexPrintf("xnow %f\n", -p[n1]);
             // mexPrintf("xnow %f",xnow);
             // The following loop calculate p_n(x) using recursive formula
      
             for(int j=1; j<n1;j++){
                 p[j+1] = c1[j]*xnow*p[j]-c2[j]*p[j-1];
             }
                
        // The derivative of p_n(x) can be calculated from p_n(x) and p_n-1(x)
            pdir = (double)n * (p[n-1]-xnow*p[n1-1])/(1-xnow*xnow);
            
            if(abs(xinc)< eps || ncount > nmax){
                found = 1;
                x[k] = xnow;
                x[i] = -xnow;
                w[k] = 2/(1-x[k]*x[k])/(pdir*pdir);
                w[i] = w[k];
            }
            xinc = -p[n1-1]/pdir;
           // mexPrintf("xinc %f\n" ,xinc);
            xnow = xnow + xinc;
         }
        }
     
   //mxFree(c1);
   //mxFree(c2);
   //mxFree(p);
}

void noThread(void);

long getNextVal(int* Vals){
long returnVal;
returnVal = CurrentIndex;
if(CurrentIndex+block<numRows){ // Will it know how many rows?
	//CurrentIndex = CurrentIndex+block-1;
    CurrentIndex = CurrentIndex + block;
	*Vals = block;
}
else{
	*Vals = numRows-CurrentIndex;
	//mexPrintf("Thread Running %d, %d, %d\n", numRows, CurrentIndex,*Vals);
	CurrentIndex = numRows;
}
return returnVal;
}

void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{ 
	double *outArray;
	DWORD dwThreadId;
	int i, nThreads = MAX_THREADS;
	HANDLE  hThreadArray[MAX_THREADS]; 

    if (nrhs != 11) {
        printf(" Value of errno: %d\n ", errno);
        return 0;
    }

	numRows = mxGetM(prhs[7]); //total rows
	CurrentIndex = 0;
	
	plhs[0] = mxCreateDoubleMatrix(numRows, 1, mxREAL); //mxReal is our data-type
    plhs[1] = mxCreateDoubleMatrix(numRows, 1, mxREAL); //mxReal is our data-type
    plhs[2] = mxCreateDoubleMatrix(numRows, 1, mxREAL); //mxReal is our data-type
    plhs[3] = mxCreateDoubleMatrix(numRows, 1, mxREAL); //mxReal is our data-type
    
	Gvv = mxGetPr(plhs[0]);	
    Gvh = mxGetPr(plhs[1]);	
    Ghv = mxGetPr(plhs[2]);	
    Ghh = mxGetPr(plhs[3]);	

    // Setup the data input matrices
	epsr = prhs[0];
    UPtr = mxGetPr(prhs[1]);
	QPtr = mxGetPr(prhs[2]);
    kPtr = mxGetPr(prhs[3]);
    mu2Ptr = mxGetPr(prhs[4]);
	mc2Ptr = mxGetPr(prhs[5]);
    kdPtr = mxGetPr(prhs[6]);
	thi = mxGetPr(prhs[7]);
    ths = mxGetPr(prhs[8]);
    phi = mxGetPr(prhs[9]);	
    phs = mxGetPr(prhs[10]);
    //psi = mxGetPr(prhs[11]);	    

	xu = mxCalloc(64, sizeof(double));
    wu = mxCalloc(64, sizeof(double));
	xc = mxCalloc(64, sizeof(double));
    wc = mxCalloc(64, sizeof(double));
    
    gauss(64,xu,wu);
    gauss(64,xc,wc);
    
    //return;
  ghMutex = CreateMutex( 
        NULL,              // default security attributes
        FALSE,             // initially not owned
        NULL);             // unnamed mutex

    nThreads = numRows/block + (numRows % block != 0);
	//Set the global count to number of threads
    if(nThreads>MAX_THREADS){
        nThreads = MAX_THREADS;    
    }
	nGlobalCount = nThreads;
	//Start the threads
    //mexPrintf("threads %d\n",nThreads);

	for (i=1; i<= nThreads; i++) {
      //  mexPrintf("creating threads\n");
		countThread[i-1]=0;
		hThreadArray[i-1] = CreateThread(NULL, //Choose default security
						 0, //Default stack size
						(LPTHREAD_START_ROUTINE)&ThreadProc,
						//Routine to execute
						(LPVOID) &i, //Thread parameter
						0, //Immediately run the thread
						&dwThreadId); //Thread Id
		
		if (hThreadArray[i-1] == NULL)
		{
			mexPrintf("Error Creating Thread#: %d\n",i);
			return;
		}
	}

    WaitForMultipleObjects(nThreads, hThreadArray, TRUE, INFINITE);
    for(i=0; i<MAX_THREADS; i++)
    {
		if (i<nThreads){
        CloseHandle(hThreadArray[i]);
		}
    } 
    
   //mxFree(xu);
   //mxFree(wu);
   //mxFree(xc);
   //mxFree(wc);
   return;
}

//Thread Routine
DWORD ThreadProc (LPVOID lpdwThreadParam ) 
{
	int index;
	int i,jj;
	int Vals=0;
	long useIndex,processIndex;
	DWORD dwWaitResult; 
    double deg = 57.295779513082320;
    double Rad = 0.017453292519943;
    double pi =  3.141592653589793;
    double two_pi = 6.283185307179586;
    double one_divided_two_pi = 0.159154943091895;
    double phir, phsr, thir, thsr, sumax;
    double ci, si, ss,cs,phsi,sphi,cphi,Ds,Di,D,cphs,sphs,coti,su,sc,phn,cpn,cpn2,spn,spn2,ctn,stn,qxy2,Ptn,hxi,hyi,hzi,y,x,aphi,xi,yi,athi,hxs,hys,hzs,aphs,xs,ys,aths,aphsi;
    double vivi,vihi,hivi,hihi,vsvs,vshs,hshs,hsvs,fvv,fvh,fhh,fhv,csi,css,sys,syi,qxy,K,Y1,Y,VV,VH,HH,HV,px,gtf,sumin,shh,svh,shv,svv,Xi,Xs,Yi,Ys,cos_two_PSIR;
    double shhR,shhI,svhR,svhI,shvR,shvI,svvR,svvI;
    double fvvR,fvvI,fvhR,fvhI,fhhR,fhhI,fhvR,fhvI;
    double U,Q,k,mu2,mc2,kd;
    bool Cx,Cy,C;
    double ngu = 64;
	double ngc = 64;
    double GvvTemp,GvhTemp,GhvTemp,GhhTemp;
   // int nu, nc;
    psi = phi;
    cos_two_PSIR = cos(2* *psi * Rad);
	
   U = *UPtr;
   Q = *QPtr;
   k = *kPtr;
   mu2 = *mu2Ptr;
   mc2 = *mc2Ptr;
   kd = *kdPtr;
    
	while(CurrentIndex<numRows){
		//reserve mutex to get the next block of data to process
		do
		{
			dwWaitResult = WaitForSingleObject( 
				ghMutex,    // handle to mutex
				INFINITE);  // no time-out interval
		}while(dwWaitResult != WAIT_OBJECT_0) ;	
		useIndex = getNextVal(&Vals); //get next block index
		//ReleaseMutex(ghMutex); //release mutex to allow for another thread to get a block
		for(jj=0;jj<Vals;jj++){       
			//countThread[index] = countThread[index]++;
			processIndex = useIndex+jj;
            
            phir=*(phi+processIndex)*Rad;
            cphi=cos(phir);
            sphi=sin(phir);
            
            phsr=*(phs+processIndex)*Rad;
            cphs=cos(phsr);
            sphs=sin(phsr);
            
            //Large scale
            thir=*(thi+processIndex)*Rad;
            ci=cos(thir);
            si=sin(thir);
            coti=1/tan(thir);
            thsr=*(ths+processIndex)*Rad;

            //Integration limits over smu2 (upwind)Eq. 29a
            sumax=6*sqrt(mu2);

            // Integration limits over smc2 (cross wind). Eq. 29b
            scmax=6*sqrt(mc2);
            scmin=-scmax;            
          
            if(coti<sumax){
                sumin=-coti;
            }
            else{
                sumin=-sumax;
            }
                
            cs=cos(thsr);
            ss=sin(thsr);
            *Gvv = 0;
            *Gvh = 0;
            *Ghv = 0;
            *Ghh = 0;
            GvvTemp =0;
            GvhTemp =0;
            GhvTemp =0;
            GhhTemp =0;
            
            for(int nu = 0;nu<64;nu++){
                su=((sumax-sumin)*xu[nu]+(sumax+sumin))/2.0;
               
                for(int nc = 0;nc<64;nc++){
                    sc=((scmax-scmin)*xc[nc]+(scmax+scmin))/2.0;
                    phn=atan2(sc,su);   // Eq. 30 - 31b
                    cpn=cos(phn);
                    cpn2=cpn*cpn;
                    spn=sin(phn);
                    spn2=spn*spn;
                    ctn=1.0/sqrt(su*su + sc*sc+1.0);
                    stn=ctn*(su*cpn + sc*spn);  

                    qxy2= 0.5 * (su * su / mu2 + sc * sc / mc2);
                    Ptn=1.0/(two_pi * sqrt(mu2 * mc2)) * exp(-qxy2);     // Slope probability distribution (41)

                    // local incident Eqs. (33a) - (33c)
                    hxi=si*sphi-ci*sc;
                    hyi=ci*su-si*cphi;
                    hzi=si*(sphi*su-cphi*sc);
                    Di=sqrt(hxi*hxi+hyi*hyi+hzi*hzi);
                    
                    y=si*sin(phir-phn);
                    x=si*ctn*cos(phir-phn)-ci*stn;
                    aphi=atan2(y,x);             // local azimuth incidence angle
                    
                    //%xi=(si*stn*cos(phir-phn)+ci*ctn);
                    xi=ctn*(si*(su*cphi+sc*sphi)+ci);
                    yi=x*cos(aphi)+y*sin(aphi);
                    
                    athi=atan2(yi,xi);   // local zenith incidence angle

                    // Local scattering  Eqs. (32a) - (32c)
                    hxs=ss*sphs+cs*sc;
                    hys=-cs*su-ss*cphs;
                    hzs=ss*(sphs*su-cphs*sc);
                    Ds=sqrt(hxs*hxs+hys*hys+hzs*hzs);

                    y=ss*sin(phsr-phn);
                    x=ss*ctn*cos(phsr-phn)+cs*stn;
                    aphs=atan2(y,x);        // local scattering azimuth angle

                    //%xs=-ss*stn*cos(phsr-phn)+cs*ctn;
                    xs=-ctn*(ss*(su*cphs+sc*sphs)-cs);
                    ys=x*cos(aphs)+y*sin(aphs);
                    aths=atan2(ys,xs);      // local scattered zenith angle

                    aphsi=aphs-aphi;

                    // Local small perturbation bistatic scattering coefficients
                    // Eqs. (38a) - (38d)                   

                    small_pert_model(&shhR,&shhI,&svhR,&svhI,&shvR,&shvI,&svvR,&svvI,epsr,athi,aths,aphi,aphs);

                    //sea_sur_small_scale(&shh,&svh,&shv,&svv,epsr,athi,aths,aphi,aphs);
                    //[shh,svh,shv,svv]=sea_sur_small_scale(epsr,athi,aths,aphi,aphs);
                    
                    // Vector cross products 
                     if(Di==0){
                       vivi=1;
                       vihi=0;
                       hivi=0;
                       hihi=1;
                     }
                     else {          // Eqs. (36a) - (36d)
                        vivi=(ctn*(-ci*(su*cphi+sc*sphi)+si)/Di);
                        vihi=(ctn*(sc*cphi-su*sphi)/Di);
                        hivi=(-(ci*(hxi*cphi+hyi*sphi)+hzi*si)/Di);
                        hihi=((hyi*cphi-hxi*sphi)/Di);
                    }
                    
                    if(Ds==0){
                        vsvs=1.0;
                        vshs=0.0;
                        hsvs=0.0;
                        hshs=1.0;
                    }
                   else{        // Eqs. (37a) - (37d)
                        vsvs=(ctn*(cs*(su*cphs+sc*sphs)+ss)/Ds);
                        vshs=((cs*(hxs*cphs+hys*sphs)-hzs*ss)/Ds);
                        hsvs=(ctn*(sc*cphs-su*sphs)/Ds);
                        hshs=((hys*cphs-hxs*sphs)/Ds);
                   }

                    //       Eqs. (39a) - (39d) 
             
                    //fvv=(vsvs*svv+vshs*shv)*vivi+(vsvs*svh+vshs*shh)*hivi;
                    fvvR =(vsvs*svvR+vshs*shvR)*vivi+(vsvs*svhR+vshs*shhR)*hivi;
                    fvvI =(vsvs*svvI+vshs*shvI)*vivi+(vsvs*svhI+vshs*shhI)*hivi;
                    
                    //fvh=(vsvs*svv+vshs*shv)*vihi+(vsvs*svh+vshs*shh)*hihi;
                    fvhR =(vsvs*svvR+vshs*shvR)*vihi+(vsvs*svhR+vshs*shhR)*hihi;
                    fvhI =(vsvs*svvI+vshs*shvI)*vihi+(vsvs*svhI+vshs*shhI)*hihi;
                    
                    //fhh=(hsvs*svv+hshs*shv)*vihi+(hsvs*svh+hshs*shh)*hihi;
                    fhhR =(hsvs*svvR+hshs*shvR)*vihi+(hsvs*svhR+hshs*shhR)*hihi;
                    fhhI =(hsvs*svvI+hshs*shvI)*vihi+(hsvs*svhI+hshs*shhI)*hihi;
                    
                    //fhv=(hsvs*svv+hshs*shv)*vivi+(hsvs*svh+hshs*shh)*hivi;                   
                    fhvR =(hsvs*svvR+hshs*shvR)*vivi+(hsvs*svhR+hshs*shhR)*hivi;
                    fhvI =(hsvs*svvI+hshs*shvI)*vivi+(hsvs*svhI+hshs*shhI)*hivi;
                    
                    // Getting Eq. 43
                    csi = cos(athi);
                    css = cos(aths);
                    csi = ctn*(si*(su*cphi+sc*sphi)+ci);
                    css = -ctn*(ss*(su*cphs+sc*sphs)-cs);
                    sys = sin(aths);
                    syi = sin(athi);

                    qxy=syi*syi+sys*sys-2.0*syi*sys*cos(aphs-aphi);
                    K= k * sqrt(qxy);        // Equation (43)
                    sea_sur_spectra(&SK,&DK,U,K,Q);

                    // Sea surface height spectrum (Annex D)
                    if(K >= kd){
                        
                        Y=(one_divided_two_pi * (1.0 + DK * cos_two_PSIR) * SK / K);   // A factor 1/K is introduced based on Eq. D-2
                        Y1=(2.0 * k * k * csi * css);
                        Y=(2.0*Y1*Y1*Y);    // Required in Eq. (42)
                        Y=(two_pi*Y)*(K>= kd);       // Based on Discussion with Prof. Johson
                    }
                    else{
                        Y=0;
                     }
                    
                    //VV=abs(fvv*fvv)*Y;     // Eqs. (42)
                    //VH=abs(fvh*fvh)*Y;
                    //HV=abs(fhv*fhv)*Y;
                    //HH=abs(fhh*fhh)*Y;
                    VV = complex_abs(fvvR*fvvR-fvvI*fvvI,fvvR*fvvI+fvvI*fvvR)*Y;
                    VH = complex_abs(fvhR*fvhR-fvhI*fvhI,fvhR*fvhI+fvhI*fvhR)*Y;
                    HV = complex_abs(fhvR*fhvR-fhvI*fhvI,fhvR*fhvI+fhvI*fhvR)*Y;
                    HH = complex_abs(fhhR*fhhR-fhhI*fhhI,fhhR*fhhI+fhhI*fhhR)*Y;

                    //   Getting Eq. (40)
                   
                    if(csi <0){
                         Xi=0;
                     }
                     else{
                        Xi=1;
                     }
                     
                    if(css <0){
                         Xs=0;
                    }
                   else{
                        Xs=1;
                   }
                    Cx=Xi&&Xs;
                    //
                     if(syi <0){
                         Yi=0;
                     }
                     else{
                        Yi=1;
                    }
                    
                    if(sys <0){
                         Ys=0;
                    }   
                    else{
                        Ys=1;
                    }
                    
                    Cy=Yi&&Ys;
                    C=Cx&&Cy;
                    
                    if(C){
                        px= (1 + si/ci*su);    // Eq. (42)
                    }
                     else{
                         px=0;
                    }

                    gtf=(sumax-sumin)*(scmax-scmin)*Ptn*wu[nu]*wc[nc]*px/4.0;

                    GvvTemp +=gtf*VV;
                    GvhTemp +=gtf*VH;
                    GhvTemp +=gtf*HV;
                    GhhTemp +=gtf*HH;                     
                }
            }
           
            // *(Gvv+processIndex) = ((sumax-sumin)*xu[0]+(sumax+sumin))/2.0;
           *(Gvv+processIndex) = GvvTemp;
           *(Gvh+processIndex) = GvhTemp;
           *(Ghv+processIndex) = GhvTemp;
           *(Ghh+processIndex) = GhhTemp;
		}
        ReleaseMutex(ghMutex); //release mutex to allow for another thread to get a block

	}
	return 0;
}

void noThread(void){
	int index;
	int i,jj;
	int Vals=0;
	long useIndex,processIndex;
	DWORD dwWaitResult; 
    double deg = 57.295779513082320;
    double A_Xval, A_Yval, A_Zval, B_Xval, B_Yval, B_Zval;
	
		for(jj=0;jj<numRows;jj++){
			/*
			
            A_Xval=*(A_X++);
            A_Yval=*(A_Y++);
            A_Zval=*(A_Z++);
            B_Xval=*(B_X++);
            B_Yval=*(B_Y++);
            B_Zval=*(B_Z++);
            *(THETA++) =  acos((A_Xval*B_Xval+A_Yval*B_Yval+A_Zval*B_Zval) / (sqrt((A_Xval*A_Xval + A_Yval*A_Yval + A_Zval*A_Zval)*(B_Xval*B_Xval + B_Yval*B_Yval + B_Zval*B_Zval)))) * deg;
		*/
		}

}