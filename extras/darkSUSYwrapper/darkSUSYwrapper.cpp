//============================================================================
// Name        : darkSUSYwrapper.cpp
// Author      : Jonathan Cornell
// Description : An exploration of how to use darkSUSY functions in C++ by
//	linking against static darkSUSY libraries at compile time. Compile this
//	with something like:
//
//	g++ darkSUSYwrapper.cpp -Wall -o "darkSUSYwrapper.o" -c
//
//  and then link against darkSUSY libraries with:
//
//  g++ -o "darkSUSYwrapper" ./src/darkSUSYwrapper.o -L $DS_PATH/lib
//		-ldarksusy -lgfortran(or whatever library is appropriate based on your
//		darkSUSY build) -lFH -lHB
//============================================================================

#include <iostream>
using namespace std;

/* The below gives the program access to the necessary common blocks and functions from the darkSUSY library */

extern "C" {void dsinit_(); void dssusy_(int *unphys, int *hwarning);
			double dsrdomega_(int *omtype, int *fast, double *xf, int *ierr, int *iwar, int *nfc);
			struct{
				double tanbe, mu, m2, m1, m3, ma;
				double mass2u[3], mass2q[3], mass2d[3], mass2l[3], mass2e[3];
				double asoftu[3], asoftd[3], asofte[3];
			} mssmpar_;
			struct{
				int prtlevel, lulog, luerr, luout;
			} dsio_;

};

int main() {

	double au3, ad3, a3half, ae3;
	int i;
	int unphys, hwarning, ierr, iwar, nfc;
	int omtype, fast;
	double xf,oh2;

	/* The following could also be pulled from darkSUSY*/
	const double mt = 172.900;
	const double mb = 4.190;
	const double me = 0.000511;
	const double mmu= .105;
	const double mtau= 1.777;
	const double vev= 246;

/*	m1=m2=m3=1000;
	mu=400;
	mA=1000;
	tanb=10;
	mQ1=mQ2=mQ3=mu1=mu2=mu3=md1=md2=md3=mL1=mL2=mL3=me1=me2=me3=2000; */
	au3=ad3=a3half=ae3=1;

	dsinit_();
	dsio_.prtlevel = 10;

	/* Here I just have hardcoded in a set of MSSM parameters that gives me a physical model */
	mssmpar_.m1=500;
	mssmpar_.m2=1000;
	mssmpar_.m3=3500;
	mssmpar_.mu=400;
	mssmpar_.ma=1000;
	mssmpar_.tanbe=10;
	for(i=0; i<=2; i++){
		mssmpar_.mass2u[i]=mssmpar_.mass2q[i]=mssmpar_.mass2d[i]=2000*2000;
		mssmpar_.mass2e[i]=mssmpar_.mass2l[i]=2000*2000;
	}

	/* The below always equal zero in MSSM 25 */
	for(i=0; i<=1; i++){
		mssmpar_.asoftu[i]=0;
		mssmpar_.asoftd[i]=0;
	}

	/* The actual MSSM 25 from arXiv:1201.0844 looks something like the commented part below.
	   Since I knew the uncommented values for the trilinear couplings gave me a physical
	   model, I used them instead */
	   
/*	mssmpar_.asofte[0] = a3half*-1*me/vev;
	mssmpar_.asofte[1] = a3half*-1*mmu/vev;
	mssmpar_.asoftu[2] = au3*-1*mt/vev;
	mssmpar_.asoftd[2] = ad3*-1*mb/vev;
	mssmpar_.asofte[2] = ae3*-1*mtau/vev; */

	mssmpar_.asofte[0] = 0;
	mssmpar_.asofte[1] = 0;
	mssmpar_.asoftu[2] = au3;
	mssmpar_.asoftd[2] = ad3;
	mssmpar_.asofte[2] = 0;

	/* Once all the above has been defined, you can just call darkSUSY functions: */
	dssusy_(&unphys, &hwarning);
	omtype=fast=1;
	oh2 = dsrdomega_(&omtype,&fast,&xf,&ierr,&iwar,&nfc);
	cout << "Omega h^2 with coannihilations is " << oh2 << endl;
	return 0;
}
