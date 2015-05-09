
namespace Gambit
{
    struct slha_leshouches1_hdec_type
    {
        Fstring<100> spinfo1;
        Fstring<100> spinfo2;
        Fstring<100> modselval;
        Fstring<20> mincom;
        Fstring<20> extcom;
        Fstring<20> softcom;
        Fstring<20> hmixcom;
    };
    struct slha_leshouches2_hdec_type
    {
        Farray< Fdouble,1,20 > minval;
        Farray< Fdouble,0,100 > extval;
        Farray< Fdouble,1,20 > smval;
        Farray< Fdouble,1,50 > massval;
        Farray< Fdouble,1,4,1,4 > nmixval;
        Farray< Fdouble,1,2,1,2 > umixval;
        Farray< Fdouble,1,2,1,2 > vmixval;
        Farray< Fdouble,1,2,1,2 > stopmixval;
        Farray< Fdouble,1,2,1,2 > sbotmixval;
        Farray< Fdouble,1,2,1,2 > staumixval;
        Farray< Fdouble,1,10 > hmixval;
        Farray< Fdouble,1,3 > gaugeval;
        Farray< Fdouble,1,100 > msoftval;
        Farray< Fdouble,1,3,1,3 > auval;
        Farray< Fdouble,1,3,1,3 > adval;
        Farray< Fdouble,1,3,1,3 > aeval;
        Farray< Fdouble,1,3,1,3 > yuval;
        Farray< Fdouble,1,3,1,3 > ydval;
        Farray< Fdouble,1,3,1,3 > yeval;
        Fdouble alphaval;
        Farray< Fdouble,1,20 > qvalue;
        Farray< Finteger,1,2 > imod;
    };
    struct widtha_hdec_type
    {
        Fdouble abrb;
        Fdouble abrl;
        Fdouble abrm;
        Fdouble abrs;
        Fdouble abrc;
        Fdouble abrt;
        Fdouble abrg;
        Fdouble abrga;
        Fdouble abrzga;
        Fdouble abrz;
        Fdouble awdth;
    };
    struct widthhl_hdec_type
    {
        Fdouble hlbrb;
        Fdouble hlbrl;
        Fdouble hlbrm;
        Fdouble hlbrs;
        Fdouble hlbrc;
        Fdouble hlbrt;
        Fdouble hlbrg;
        Fdouble hlbrga;
        Fdouble hlbrzga;
        Fdouble hlbrw;
        Fdouble hlbrz;
        Fdouble hlbra;
        Fdouble hlbraz;
        Fdouble hlbrhw;
        Fdouble hlwdth;
    };
    struct widthhh_hdec_type
    {
        Fdouble hhbrb;
        Fdouble hhbrl;
        Fdouble hhbrm;
        Fdouble hhbrs;
        Fdouble hhbrc;
        Fdouble hhbrt;
        Fdouble hhbrg;
        Fdouble hhbrga;
        Fdouble hhbrzga;
        Fdouble hhbrw;
        Fdouble hhbrz;
        Fdouble hhbrh;
        Fdouble hhbra;
        Fdouble hhbraz;
        Fdouble hhbrhw;
        Fdouble hhwdth;
    };
    struct widthhc_hdec_type
    {
        Fdouble hcbrb;
        Fdouble hcbrl;
        Fdouble hcbrm;
        Fdouble hcbrbu;
        Fdouble hcbrs;
        Fdouble hcbrc;
        Fdouble hcbrt;
        Fdouble hcbrw;
        Fdouble hcbra;
        Fdouble hcwdth;
    };
    struct wisusy_hdec_type
    {
        Farray< Fdouble,1,2,1,2 > hlbrsc;
        Farray< Fdouble,1,4,1,4 > hlbrsn;
        Farray< Fdouble,1,2,1,2 > hhbrsc;
        Farray< Fdouble,1,4,1,4 > hhbrsn;
        Farray< Fdouble,1,2,1,2 > habrsc;
        Farray< Fdouble,1,4,1,4 > habrsn;
        Farray< Fdouble,1,2,1,4 > hcbrsu;
        Fdouble hlbrcht;
        Fdouble hhbrcht;
        Fdouble habrcht;
        Fdouble hlbrnet;
        Fdouble hhbrnet;
        Fdouble habrnet;
        Fdouble hcbrcnt;
        Fdouble hlbrsl;
        Fdouble hhbrsl;
        Fdouble hcbrsl;
        Fdouble habrsl;
        Fdouble habrst;
        Fdouble habrsb;
        Fdouble hhbrsq;
        Farray< Fdouble,1,2,1,2 > hhbrst;
        Farray< Fdouble,1,2,1,2 > hhbrsb;
        Fdouble hhbrsqt;
        Fdouble hcbrsq;
        Farray< Fdouble,1,2,1,2 > hcbrstb;
        Fdouble hcbrsqt;
        Fdouble hlbrsq;
        Fdouble hlbrsqt;
    };
    struct wisfer_hdec_type
    {
        Fdouble bhlslnl;
        Fdouble bhlslel;
        Fdouble bhlsler;
        Fdouble bhlsqul;
        Fdouble bhlsqur;
        Fdouble bhlsqdl;
        Fdouble bhlsqdr;
        Fdouble bhlst;
        Fdouble bhlsb;
        Fdouble bhlstau;
        Fdouble bhhslnl;
        Fdouble bhhslel;
        Fdouble bhhsler;
        Fdouble bhhsqul;
        Fdouble bhhsqur;
        Fdouble bhhsqdl;
        Fdouble bhhsqdr;
        Fdouble bhhst;
        Fdouble bhhsb;
        Fdouble bhhstau;
        Fdouble bhastau;
        Fdouble bhasb;
        Fdouble bhast;
        Fdouble bhcsl00;
        Fdouble bhcsl11;
        Fdouble bhcsl21;
        Fdouble bhcsq;
        Fdouble bhcstb;
    };
    struct susyhitin_type
    {
        Fdouble flagshsin;
        Fdouble amsin;
        Fdouble amcin;
        Fdouble ammuonin;
        Fdouble alphin;
        Fdouble gamwin;
        Fdouble gamzin;
        Fdouble vusin;
        Fdouble vcbin;
        Fdouble rvubin;
    };
}
