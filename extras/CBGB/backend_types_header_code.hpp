
namespace Gambit
{
    struct pacodes_type
    {
        Farray< Finteger,1,2 > kse;
        Farray< Finteger,1,2 > ksmu;
        Farray< Finteger,1,2 > kstau;
        Farray< Finteger,1,2 > ksu;
        Farray< Finteger,1,2 > ksd;
        Farray< Finteger,1,2 > ksc;
        Farray< Finteger,1,2 > kss;
        Farray< Finteger,1,2 > kst;
        Farray< Finteger,1,2 > ksb;
        Farray< Finteger,1,4 > kn;
        Farray< Finteger,1,2 > kcha;
        Farray< Finteger,1,3 > knu;
        Farray< Finteger,1,3 > kl;
        Farray< Finteger,1,3 > kqu;
        Farray< Finteger,1,3 > kqd;
        Farray< Finteger,1,6 > ksnu;
        Farray< Finteger,1,6 > ksl;
        Farray< Finteger,1,6 > ksqu;
        Farray< Finteger,1,6 > ksqd;
        Fstring<8> pname;
    };
    struct mspctm_type
    {
        Farray< Freal8,0,50 > mass;
        Farray< Freal8,0,50 > runmass;
        Freal8 mu2gev;
        Freal8 md2gev;
        Freal8 ms2gev;
        Freal8 mcmc;
        Freal8 mbmb;
        Freal8 mtmt;
    };
    struct widths_type
    {
        Farray< Freal8,0,50 > width;
    };
    struct intdof_type
    {
        Farray< Finteger,0,50 > kdof;
    };
    struct vrtxs_type
    {
        Farray< Fcomplex16,1,50,1,50,1,50 > gl;
        Farray< Fcomplex16,1,50,1,50,1,50 > gr;
    };
    struct smruseful_type
    {
        Freal8 s2thw;
        Freal8 sinthw;
        Freal8 costhw;
        Freal8 delrho;
        Freal8 alph3mz;
        Freal8 gfermi;
        Freal8 s2wmz;
        Freal8 swmz;
        Freal8 cwmz;
    };
    struct smcuseful_type
    {
        Fstring<5> roption;
    };
    struct couplingconstants_type
    {
        Freal8 g2weak;
        Freal8 gyweak;
        Freal8 g3stro;
        Freal8 alphem;
        Freal8 alph3;
        Farray< Freal8,1,12 > yukawa;
        Freal8 g2wmz;
        Freal8 gywmz;
    };
    struct sckm_type
    {
        Freal8 ckms12;
        Freal8 ckms23;
        Freal8 ckms13;
        Freal8 ckmdelta;
    };
    struct mixing_type
    {
        Farray< Fcomplex16,1,3,1,3 > ckm;
    };
    struct mssmtype_type
    {
        Finteger modeltype;
    };
    struct mssmpar_type
    {
        Freal8 tanbe;
        Freal8 mu;
        Freal8 m2;
        Freal8 m1;
        Freal8 m3;
        Freal8 ma;
        Farray< Freal8,1,3 > mass2u;
        Farray< Freal8,1,3 > mass2q;
        Farray< Freal8,1,3 > mass2d;
        Farray< Freal8,1,3 > mass2l;
        Farray< Freal8,1,3 > mass2e;
        Farray< Freal8,1,3 > asoftu;
        Farray< Freal8,1,3 > asoftd;
        Farray< Freal8,1,3 > asofte;
    };
    struct mssmswitch_type
    {
        Freal8 msquarks;
        Freal8 msleptons;
        Finteger higloop;
        Finteger neuloop;
        Finteger bsgqcd;
        Finteger higwid;
    };
    struct sfermionmass_type
    {
        Farray< Freal8,1,3 > massup1;
        Farray< Freal8,1,3 > massup2;
        Farray< Freal8,1,3 > thetamixu;
        Farray< Freal8,1,3 > massdn1;
        Farray< Freal8,1,3 > massdn2;
        Farray< Freal8,1,3 > thetamixd;
        Farray< Freal8,1,3 > masssn;
        Farray< Freal8,1,3 > masssl1;
        Farray< Freal8,1,3 > masssl2;
        Farray< Freal8,1,3 > thetamixsl;
    };
    struct mssmwidths_type
    {
        Farray< Freal8,1,32,1,4 > hdwidth;
    };
    struct mssmmixing_type
    {
        Farray< Fcomplex16,1,4,1,4 > neunmx;
        Farray< Fcomplex16,1,2,1,2 > chaumx;
        Farray< Fcomplex16,1,2,1,2 > chavmx;
        Farray< Fcomplex16,1,3,1,3 > slulmx;
        Farray< Fcomplex16,1,6,1,3 > sldlmx;
        Farray< Fcomplex16,1,6,1,3 > sldrmx;
        Farray< Fcomplex16,1,6,1,3 > squlmx;
        Farray< Fcomplex16,1,6,1,3 > squrmx;
        Farray< Fcomplex16,1,6,1,3 > sqdlmx;
        Farray< Fcomplex16,1,6,1,3 > sqdrmx;
        Freal8 alpha;
        Freal8 mix_stop;
        Freal8 mix_sbot;
        Freal8 mix_stau;
    };
}
