/************************************************************************//**
 * File: srwlclient.cpp
 * Description: Demo C/C++ client 
 * Project: Synchrotron Radiation Workshop Library (SRWLib)
 * First release: October 2010
 *
 * SRW is Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
 * SRW C/C++ API (SRWLIB) is Copyright (C) European XFEL, Hamburg, Germany
 * All Rights Reserved
 *
 * @author O.Chubar, G.Geloni, L.Samoylova
 * @version 0.04
 ***************************************************************************/

#include "srwlib.h"
#include <iostream>
#include <fstream>
#include <cstring> //necessary for strcpy, etc...
#include <cstdlib> //necessary for atoi
#include <sstream>
#include <iomanip>
#include <random>
#include <chrono>
#include <mcpl.h>

using namespace std;

/************************************************************************//**
 * Auxiliary function dedicated to process errors reported by Library
 ***************************************************************************/
void ProcRes(int er)
{
	char ErrorBuf[2048];

	if(er == 0) return;
	else
	{
		srwlUtiGetErrText(ErrorBuf, er);

		cout << endl;
		if(er < 0) 
		{//Print Warning:
			cout << "WARNING: " << ErrorBuf << endl;
		}
		else 
		{//Just print Error Message:
			cout << "ERROR: " << ErrorBuf << endl;
			//cout << "Press Enter to exit" << endl;
			//getchar();
			//exit(0);
		}
	}
}

/************************************************************************//**
 * Wavefront modification (re-allocation) function; to be called by pointer from SRWLIB
 ***************************************************************************/
int ModifySRWLWfr(int action, SRWLWfr* pWfr, char pol)
{
	if(pWfr == 0) return -1; //returning non-zero means Wfr modification did not succeed; no throwing allowed here
	if((action < 0) || (action > 2)) return -1;

	long numTot = pWfr->mesh.ne*pWfr->mesh.nx*pWfr->mesh.ny*2;
	if(numTot <= 0) return 0; //or delete the previous array, still?

	int ExNeeded = ((pol == 0) || (pol == 'x') || (pol == 'X'))? 1 : 0;
	int EyNeeded = ((pol == 0) || (pol == 'y') || (pol == 'Y') || (pol == 'z') || (pol == 'Z'))? 1 : 0;

		//DEBUG
		cout << "nx= " << pWfr->mesh.nx << ", ny= " << pWfr->mesh.ny << endl; 

	if(action == 0) 
	{//just delete existing wavefront data
		if(ExNeeded)
		{
			if(pWfr->arEx) delete[] pWfr->arEx;
			pWfr->arEx = 0;
			if(pWfr->arMomX) delete[] pWfr->arMomX;
			pWfr->arMomX = 0;
		}
		if(EyNeeded)
		{
			if(pWfr->arEy) delete[] pWfr->arEy;
			pWfr->arEy = 0;
			if(pWfr->arMomY) delete[] pWfr->arMomY;
			pWfr->arMomY = 0;
		}
	}
	else if(action == 1)
	{//allocate new wavefront data (without checking/deleting any existing data)
		if(ExNeeded)
		{
			pWfr->arEx = (char*)(new float[numTot]);
			pWfr->arMomX = new double[11*pWfr->mesh.ne];
		}
		if(EyNeeded)
		{
			pWfr->arEy = (char*)(new float[numTot]);
			pWfr->arMomY = new double[11*pWfr->mesh.ne];
		}
	}
	else if(action == 2)
	{//modify wavefront size (numbers of points vs photon energy, horizontal or vertical position)
		if(ExNeeded)
		{//using realloc could perhaps be more efficient here
			if(pWfr->arEx) delete[] pWfr->arEx;
			pWfr->arEx = (char*)(new float[numTot]);
			if(pWfr->arMomX) delete[] pWfr->arMomX;
			pWfr->arMomX = new double[11*pWfr->mesh.ne];
		}
		if(EyNeeded)
		{//using realloc could perhaps be more efficient here
			if(pWfr->arEy) delete[] pWfr->arEy;
			pWfr->arEy = (char*)(new float[numTot]);
			if(pWfr->arMomY) delete[] pWfr->arMomY;
			pWfr->arMomY = new double[11*pWfr->mesh.ne];
		}
	}
	return 0;
}

/************************************************************************//**
 * Auxiliary function to read line from file and to extract a number from it
 ***************************************************************************/
template<class T> void AuxReadLineAndExtractNumber(T& numOut, ifstream& fIn, const string& sCom) //throw(...) 
{
	string sIn;
	getline(fIn, sIn); if(!fIn.good()) throw -1;

	size_t posComStart = sIn.find(sCom);
	if(posComStart != 0) throw -1;

	posComStart = sCom.length();
	size_t posComEnd = sIn.find(sCom, posComStart);

	istringstream is(sIn.substr(posComStart, posComEnd - 1));
	is >> numOut;
}

/************************************************************************//**
 * Auxiliary function to read 3D magnetic field data from ASCII file
 * File format is not flexible!
 ***************************************************************************/
int AuxReadInMagFld3D(SRWLMagFld3D* pFld, const char* strFileName)
{
	if((strFileName == 0) || (pFld == 0)) return -1;
		//cout << strFileName << endl;
	ifstream f(strFileName);
	if(!f.is_open()) return -1;

	string sRead;
	double xStart = 0, xStep = 0, yStart = 0, yStep = 0, zStart = 0, zStep = 0;
	int xNp = 1, yNp = 1, zNp = 1;
	pFld->arBx = 0; pFld->arBy = 0; pFld->arBx = 0;

	try
	{
		getline(f, sRead); if(!f.good()) return -1; //1st line: just pass 
		AuxReadLineAndExtractNumber(xStart, f, "#"); //2nd line: initial X position [m]; it will not actually be used
		AuxReadLineAndExtractNumber(xStep, f, "#"); //3rd line: step vs X [m]
		AuxReadLineAndExtractNumber(xNp, f, "#"); //4th line: number of points vs X

		AuxReadLineAndExtractNumber(yStart, f, "#"); //5th line: initial Y position [m]; it will not actually be used
		AuxReadLineAndExtractNumber(yStep, f, "#"); //6th line: step vs Y [m]
		AuxReadLineAndExtractNumber(yNp, f, "#"); //7th line: number of points vs Y

		AuxReadLineAndExtractNumber(zStart, f, "#"); //8th line: initial Z position [m]; it will not actually be used
		AuxReadLineAndExtractNumber(zStep, f, "#"); //9th line: step vs Z [m]
		AuxReadLineAndExtractNumber(zNp, f, "#"); //10th line: number of points vs Z

		pFld->rx = xStep*(xNp - 1);
		pFld->nx = xNp;
		pFld->ry = yStep*(yNp - 1);
		pFld->ny = yNp;
		pFld->rz = zStep*(zNp - 1);
		pFld->nz = zNp;

		long bNp = xNp*yNp*zNp;
		pFld->arBx = new double[bNp];
		pFld->arBy = new double[bNp];
		pFld->arBz = new double[bNp];

		double *t_arBx = pFld->arBx, *t_arBy = pFld->arBy, *t_arBz = pFld->arBz;
		for(long i=0; i<bNp; i++)
		{
			//lines from 11th: Magnetic Field Components Bx, By, Bz
			getline(f, sRead); if(!f.good()) return -1;
			istringstream is(sRead);
			is >> *(t_arBx++); 
			is >> *(t_arBy++); 
			is >> *(t_arBz++); 
		}
		f.close();
	}
	catch(int erNo) 
	{
		if(pFld->arBx != 0) { delete[] pFld->arBx; pFld->arBx = 0;}
		if(pFld->arBy != 0) { delete[] pFld->arBy; pFld->arBy = 0;}
		if(pFld->arBz != 0) { delete[] pFld->arBz; pFld->arBz = 0;}
		if(f.is_open()) f.close();
		return erNo;
	}
	return 0;
}

/************************************************************************//**
 * Auxiliary function to save trajectory data to ASCII file
 * File format is not flexible!
 ***************************************************************************/
int AuxSaveTrajData(SRWLPrtTrj* pTrj, const char* strFileName)
{
	if((strFileName == 0) || (pTrj == 0)) return -1;

	ofstream f(strFileName);
	if(!f.is_open()) return -1;

	int nCh = 12; //number of characters in each output value
	f.precision(nCh);

	f << "#ct [m], X [m], BetaX [rad], Y [m], BetaY [rad], Z [m], BetaZ [m]" << endl;

	double ctStep = (pTrj->np > 0)? (pTrj->ctEnd - pTrj->ctStart)/(pTrj->np - 1) : 0;
	double ct = pTrj->ctStart;
	double *t_arX = pTrj->arX, *t_arXp = pTrj->arXp;
	double *t_arY = pTrj->arY, *t_arYp = pTrj->arYp;
	double *t_arZ = pTrj->arZ, *t_arZp = pTrj->arZp;
	for(long i=0; i<pTrj->np; i++)
	{
		f << ct << '\t' << *(t_arX++) << '\t' << *(t_arXp++) << '\t' << *(t_arY++) << '\t' << *(t_arYp++) << '\t' << *(t_arZ++) << '\t' << *(t_arZp++) << endl;
		ct += ctStep;
	}
	f << ends;

	if(f.is_open()) f.close();
	return 0;
}

/************************************************************************//**
 * Auxiliary function to save intensity data to ASCII file
 * File format is not flexible!
 ***************************************************************************/
int AuxSaveIntensData(float* arI, double eSt, double eFi, int ne, double xSt, double xFi, int nx, double ySt, double yFi, int ny, const char* strFileName)
{
	if((strFileName == 0) || (arI == 0)) return -1;

	ofstream f(strFileName);
	if(!f.is_open()) return -1;

	int nCh = 8; //number of characters in each output value
	f.precision(nCh);
	f << "C-aligned Intensity (inner loop is vs photon energy, outer loop vs vertical position)" << endl;
	f << '#' << eSt << " #Initial Photon Energy [eV]\n";
    f << '#' << eFi << " #Final Photon Energy [eV]\n";
    f << '#' << ne << " #Number of points vs Photon Energy\n";
    f << '#' << xSt << " #Initial Horizontal Position [m]\n";
    f << '#' << xFi << " #Final Horizontal Position [m]\n";
    f << '#' << nx << " #Number of points vs Horizontal Position\n";
    f << '#' << ySt << " #Initial Vertical Position [m]\n";
    f << '#' << yFi << " #Final Vertical Position [m]\n";
    f << '#' << ny << " #Number of points vs Vertical Position\n";

	float *t_atI = arI;
	for(long i=0; i<(ne*nx*ny); i++) f << " " << *(t_atI++) << endl;
	f << ends;
	if(f.is_open()) f.close();
	return 0;
}

int gen_mcpl_particle(mcpl_particle_t *a, double p, double x, double y, double z, double e, double t, double x0,double y0,double z0){
  double dist = z-z0;
  double ux,uy,uz;
  ux=(x-x0)/dist;
  uy=(y-y0)/dist;
  uz=1.0;
  double uu = sqrt(ux*ux+uy*uy+uz*uz);
  ux/=uu;
  uy/=uu;
  uz/=uu;
  
  /*SRW energy is in eV -> MCPL in MeV*/
  a->ekin=e/1e6;

  /*MCPL wants cm - not m*/
  a->position[0]=x0*100; a->position[1]=y0*100; a->position[2]=z0*100; 
  a->direction[0]=ux; a->direction[1]=uy; a->direction[2]=uz;
  a->time=t*1e3;
  a->weight=p;
  return 0;
}


/************************************************************************//**
 * Example#3: Calculating synchrotron (undulator) radiation emitted by an electron travelling in ellipsoidal undulator
 ***************************************************************************/
int SRWLIB2MCPL(long long ncount, const char *filename)
{
	cout << "SRWLIB C Client Example:" << endl; 
	cout << "Calculating synchrotron (undulator) radiation emitted by an electron travelling in linear undulator std_sp8" << endl; 
	
	//***********Undulator
	int numPer = 140; //Number of ID Periods (without counting for terminations
	double undPer = 0.032; //Period Length [m]
	double Bx = 0; //Peak Horizontal field [T]
	double By = 0.345115; //Peak Vertical field [T]
	double phBx = 0; //Initial Phase of the Horizontal field component
	double phBy = 0; //Initial Phase of the Vertical field component
	int sBx = -1; //Symmetry of the Horizontal field component vs Longitudinal position
	int sBy = 1; //Symmetry of the Vertical field component vs Longitudinal position
	double xcID = 0; //Transverse Coordinates of Undulator Center [m]
	double ycID = 0;
	double zcID = 0; //Longitudinal Coordinate of Undulator Center [m]

	SRWLMagFldH harm;
	harm.n = 1; //harmonic number
	harm.h_or_v = 'v'; //magnetic field plane: horzontal ('h') or vertical ('v')
	harm.B = By; //magnetic field amplitude [T]
	harm.ph = 0; //phase [rad]
	harm.s = sBy; //symmetry vs longitudinal position: 1 - symmetric (B ~ cos(2*Pi*n*z/per + ph)) , -1 - anti-symmetric (B ~ sin(2*Pi*n*z/per + ph))
	harm.a = 1; //coefficient for transverse depenednce: B*cosh(2*Pi*n*a*y/per)*cos(2*Pi*n*z/per + ph)

	SRWLMagFldU und; //Ellipsoidal Undulator
	und.arHarm = &harm; //arH; //arHarmonics; //array of field harmonics
	und.nHarm = 1; //number of field harmonics
	und.per = undPer; //period length [m]
	und.nPer = numPer; //number of periods

	SRWLMagFldC magFldCnt;
	void *vArMagFld[] = {(void*)(&und)};
	magFldCnt.arMagFld = vArMagFld; //array of pointers to magnetic field elements
	magFldCnt.arMagFldTypes = "u"; //types of magnetic field elements in arMagFld array
	double auxArXcID[] = {xcID};
	double auxArYcID[] = {ycID};
	double auxArZcID[] = {zcID};
	magFldCnt.arXc = auxArXcID; //horizontal center positions of magnetic field elements in arMagFld array
	magFldCnt.arYc = auxArYcID; //vertical center positions of magnetic field elements in arMagFld array
	magFldCnt.arZc = auxArZcID; //longitudinal center positions of magnetic field elements in arMagFld array
	magFldCnt.nElem = 1; //number of magnetic field elements in arMagFld array

	//***********Electron Beam
	SRWLPartBeam elecBeam;
	elecBeam.Iavg = 0.1; //Average Current [A]
	elecBeam.partStatMom1.x = 0.; //Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
	elecBeam.partStatMom1.y = 0.;
	elecBeam.partStatMom1.z = -0.5*undPer*(numPer + 4); //Initial Longitudinal Coordinate (set before the ID)
	elecBeam.partStatMom1.xp = 0; //Initial Relative Transverse Velocities
	elecBeam.partStatMom1.yp = 0;
	elecBeam.partStatMom1.gamma = 8./0.51099890221e-03; //Relative Energy
	elecBeam.partStatMom1.relE0 = 1; //Rest mass (energy) in units of electron rest mass: =1 for electron, =1836.1526988 (=938.272013/0.510998902) for proton
	elecBeam.partStatMom1.nq = -1; //Charge of the particle related to absolute value of electron charge: =-1 for electron, =1 for positron and for proton

	//***********Precision
	int meth = 1; //SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
	double relPrec = 0.01; //relative precision
	double zStartInteg = 0; //longitudinal position to start integration (effective if < zEndInteg)
	double zEndInteg = 0; //longitudinal position to finish integration (effective if > zStartInteg)
	int npTraj = 20000;
	double sampFactNxNyForProp = 0; //sampling factor for adjusting nx, ny (effective if > 0)
	double arPrecPar[] = {meth, relPrec, zStartInteg, zEndInteg, npTraj, 0, sampFactNxNyForProp};

	//***********Wavefronts
	SRWLWfr wfr1; //For spectrum vs photon energy
	wfr1.mesh.ne = 1; //Numbers of points vs Photon Energy, Horizontal and Vertical Positions
	wfr1.mesh.nx = wfr1.mesh.ny = 1;
	wfr1.mesh.zStart = 20.; //Longitudinal Position [m] at which SR has to be calculated
	wfr1.mesh.eStart = 8000.; //Initial Photon Energy [eV]
	wfr1.mesh.eFin = 32000.; //Final Photon Energy [eV]
	wfr1.mesh.xStart = -1.0e-4; //Initial Horizontal Position [m]
	wfr1.mesh.xFin = 1.0e-4; //Final Horizontal Position [m]
	wfr1.mesh.yStart = -1.0e-5; //Initial Vertical Position [m]
	wfr1.mesh.yFin = 1.0e-5; //Final Vertical Position [m]
	wfr1.partBeam = elecBeam;
	wfr1.presCA = 0; //presentation/domain: 0- coordinates, 1- angles
	wfr1.presFT = 0; //presentation/domain: 0- frequency (photon energy), 1- time
	long numTot = wfr1.mesh.ne*wfr1.mesh.nx*wfr1.mesh.ny*2;
	float *arEx1 = new float[numTot];
	float *arEy1 = new float[numTot];
	wfr1.arEx = (char*)(arEx1); //horizontal and vertical electric field component arrays
	wfr1.arEy = (char*)(arEy1);
	wfr1.arElecPropMatr = new double[20];
	wfr1.arWfrAuxData = new double[30];
	wfr1.arMomX = new double[11*wfr1.mesh.ne];
	wfr1.arMomY = new double[11*wfr1.mesh.ne];

        SRWLWfr wfr = wfr1;
        
        //*************set up MCPL file
        mcpl_outfile_t outputfile;
        mcpl_particle_t *particle,Particle;
        int userflagenabled;
        
        char extension[128]="";
        char myfilename[2048],*stripext;
        /*strip extension from filename*/
        

        strncpy(myfilename,filename,2047);
        if ((stripext=strrchr(myfilename,'.'))){
          fprintf(stdout, "WARNING: (SRW2MCPL): the \"%s\" file extension will be overwritten by \".mcpl\"\n",stripext);
          *stripext='\0';
        }
        sprintf(extension,"mcpl");

        char line[256];
        outputfile = mcpl_create_outfile(myfilename);
        /*reset filename to be whatever mcpl actually calls it. It may have added .mcpl*/
        snprintf(myfilename,strlen(myfilename)+5,"%s",mcpl_outfile_filename(outputfile));

        snprintf(line,255,"SRW2MCPL");
        mcpl_hdr_set_srcname(outputfile,line);
        mcpl_enable_universal_pdgcode(outputfile,22);/*all particles are photons*/
        snprintf(line,255,"Output by SRW2MCPL");
        mcpl_hdr_add_comment(outputfile,line);
        /*the mcpl-file should now be set up*/


        /*pointer to the single particle storage area*/
        particle=&Particle;

        //set up random numbers here.
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();     
        std::default_random_engine generator(seed);
        std::uniform_real_distribution<double> distribution(0.0,1.0);
  
        /*Here's the particle loop*/
	float *arI1 = new float[1];
        for (int idx=0; idx<ncount; idx++){
          if(wfr1.mesh.xStart!=0){
            double xrg=wfr1.mesh.xFin-wfr1.mesh.xStart;
            wfr.mesh.xFin=wfr.mesh.xStart=distribution(generator)*xrg + wfr1.mesh.xStart;
          }
          if(wfr1.mesh.yStart!=0){
            double yrg=wfr1.mesh.yFin-wfr1.mesh.yStart;
            wfr.mesh.yFin=wfr.mesh.yStart=distribution(generator)*yrg + wfr1.mesh.yStart;
          }
          double erg=wfr1.mesh.eFin-wfr1.mesh.eStart;
          wfr.mesh.eStart=distribution(generator)*erg + wfr1.mesh.eStart;
          wfr.mesh.eFin=wfr.mesh.eStart+1.0;
          ProcRes(srwlCalcElecFieldSR(&wfr, 0, &magFldCnt, arPrecPar));
	  ProcRes(srwlCalcIntFromElecField((char*)arI1, &wfr, 6, 0, 0, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart));
         
          /*fill an MCPL particle with values the 0,0,0 starting point could be fixed*/
          double time=0;
          gen_mcpl_particle(particle,arI1[0], wfr.mesh.xStart, wfr.mesh.yStart, 20, wfr.mesh.eStart, time, 0,0,0); 
          
          mcpl_add_particle(outputfile,particle);
           
          fprintf(stdout,"%g %g\n",wfr.mesh.eStart,arI1[0]);
        }

	//**********************Deallocating memory
	delete[] arEx1;
	delete[] arEy1;
	delete[] wfr1.arElecPropMatr;
	delete[] wfr1.arWfrAuxData;
	delete[] wfr1.arMomX;
	delete[] wfr1.arMomY;
	delete[] arI1;
	  
        mcpl_closeandgzip_outfile(outputfile); 
  
        return 0;
}


/************************************************************************//**
 * Main: illustrates and tests the basic functionality of SRW Library
 ***************************************************************************/
int main(int argc, char* argv[])
{
        char ofn[512];
        long long N;
        strncpy(ofn,"voutput.mcpl",512);
        N=10000;
        switch (argc){
          case 3:
            strncpy(ofn,argv[2],512);
          case 2:
            N=strtoll(argv[1],NULL,10);
        };


	char sVersNoSRW[1024], sVersNoSRWLIB[1024];
	ProcRes(srwlUtiVerNo(sVersNoSRW, 1));
	ProcRes(srwlUtiVerNo(sVersNoSRWLIB, 2));

	cout << "SRW Version: " << sVersNoSRW << endl; 
	cout << "SRWLIB Version: " << sVersNoSRWLIB << endl; 

	//Setting "callback" function pointer
	srwlUtiSetWfrModifFunc(&ModifySRWLWfr);


	if(SRWLIB2MCPL(N, ofn)) cout << "SRW2MCPL was not executed correctly" << endl;
	return 0;
}

/***************************************************************************/
