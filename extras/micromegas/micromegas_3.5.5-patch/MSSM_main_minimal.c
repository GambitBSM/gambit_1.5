#define CLEAN    to clean intermediate files

#include"../sources/micromegas.h"
#include"../sources/micromegas_aux.h"
#include"lib/pmodel.h"

int main(int argc,char** argv)
{
    //ForceUG=1;  /* to Force Unitary Gauge assign 1 */

    printf("Initial file  \"%s\"\n",argv[1]);
    lesHinput(argv[1]);

    double pA0[2],pA5[2],nA0[2],nA5[2];
    double Nmass=0.939; /*nucleon mass*/
    double SCcoeff;        
    char cdmName[10];

    //sortOddParticles(cdmName);

    //  qNumbers(cdmName,&spin2, &charge3, &cdim);

    printf("\n==== Calculation of CDM-nucleons amplitudes  =====\n");   

    nucleonAmplitudes(FeScLoop, pA0,pA5,nA0,nA5);
    printf("CDM-nucleon micrOMEGAs amplitudes:\n");
    printf("proton:  SI  %.3E  SD  %.3E\n",pA0[0]*2,pA5[0]*2);
    printf("neutron: SI  %.3E  SD  %.3E\n",nA0[0]*2,nA5[0]*2); 

    SCcoeff=4/M_PI*3.8937966E8*pow(Nmass*Mcdm/(Nmass+ Mcdm),2.);
    printf("\n==== CDM-nucleon cross sections[pb] ====\n");
    printf(" proton  SI %.3E  SD %.3E\n",SCcoeff*pA0[0]*pA0[0],3*SCcoeff*pA5[0]*pA5[0]);
    printf(" neutron SI %.3E  SD %.3E\n",SCcoeff*nA0[0]*nA0[0],3*SCcoeff*nA5[0]*nA5[0]);

    return 0;
}
