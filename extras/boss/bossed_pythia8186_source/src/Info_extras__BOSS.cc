#include <map>
#include <string>
#include <vector>
#include <ostream>
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"
#include "Pythia8/Info.h"

void Pythia8::Info::list__BOSS() const
{
    list();
}


bool Pythia8::Info::hasSub__BOSS() const
{
    return hasSub();
}


std::basic_string<char,std::char_traits<char>,std::allocator<char> > Pythia8::Info::nameSub__BOSS() const
{
    return nameSub();
}


int Pythia8::Info::codeSub__BOSS() const
{
    return codeSub();
}


int Pythia8::Info::nFinalSub__BOSS() const
{
    return nFinalSub();
}


int Pythia8::Info::id1__BOSS() const
{
    return id1();
}


int Pythia8::Info::id2__BOSS() const
{
    return id2();
}


double Pythia8::Info::x1__BOSS() const
{
    return x1();
}


double Pythia8::Info::x2__BOSS() const
{
    return x2();
}


double Pythia8::Info::y__BOSS() const
{
    return y();
}


double Pythia8::Info::tau__BOSS() const
{
    return tau();
}


int Pythia8::Info::id1pdf__BOSS() const
{
    return id1pdf();
}


int Pythia8::Info::id2pdf__BOSS() const
{
    return id2pdf();
}


double Pythia8::Info::x1pdf__BOSS() const
{
    return x1pdf();
}


double Pythia8::Info::x2pdf__BOSS() const
{
    return x2pdf();
}


double Pythia8::Info::pdf1__BOSS() const
{
    return pdf1();
}


double Pythia8::Info::pdf2__BOSS() const
{
    return pdf2();
}


double Pythia8::Info::QFac__BOSS() const
{
    return QFac();
}


double Pythia8::Info::Q2Fac__BOSS() const
{
    return Q2Fac();
}


double Pythia8::Info::alphaS__BOSS() const
{
    return alphaS();
}


double Pythia8::Info::alphaEM__BOSS() const
{
    return alphaEM();
}


double Pythia8::Info::QRen__BOSS() const
{
    return QRen();
}


double Pythia8::Info::Q2Ren__BOSS() const
{
    return Q2Ren();
}


double Pythia8::Info::scalup__BOSS() const
{
    return scalup();
}


double Pythia8::Info::mHat__BOSS() const
{
    return mHat();
}


double Pythia8::Info::sHat__BOSS() const
{
    return sHat();
}


double Pythia8::Info::tHat__BOSS() const
{
    return tHat();
}


double Pythia8::Info::uHat__BOSS() const
{
    return uHat();
}


double Pythia8::Info::pTHat__BOSS() const
{
    return pTHat();
}


double Pythia8::Info::pT2Hat__BOSS() const
{
    return pT2Hat();
}


double Pythia8::Info::m3Hat__BOSS() const
{
    return m3Hat();
}


double Pythia8::Info::m4Hat__BOSS() const
{
    return m4Hat();
}


double Pythia8::Info::thetaHat__BOSS() const
{
    return thetaHat();
}


double Pythia8::Info::phiHat__BOSS() const
{
    return phiHat();
}


std::basic_string<char,std::char_traits<char>,std::allocator<char> > Pythia8::Info::nameProc__BOSS()
{
    return nameProc();
}


long int Pythia8::Info::nTried__BOSS()
{
    return nTried();
}


long int Pythia8::Info::nSelected__BOSS()
{
    return nSelected();
}


long int Pythia8::Info::nAccepted__BOSS()
{
    return nAccepted();
}


double Pythia8::Info::sigmaGen__BOSS()
{
    return sigmaGen();
}


double Pythia8::Info::sigmaErr__BOSS()
{
    return sigmaErr();
}


void Pythia8::Info::setCounter__BOSS(int i)
{
    setCounter(i);
}


void Pythia8::Info::addCounter__BOSS(int i)
{
    addCounter(i);
}


void Pythia8::Info::errorMsg__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > messageIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > extraIn, bool showAlways)
{
    errorMsg(messageIn, extraIn, showAlways);
}


void Pythia8::Info::errorMsg__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > messageIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > extraIn)
{
    errorMsg(messageIn, extraIn);
}


void Pythia8::Info::errorMsg__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > messageIn)
{
    errorMsg(messageIn);
}


void Pythia8::Info::errorStatistics__BOSS()
{
    errorStatistics();
}


void Pythia8::Info::setType__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int codeIn, int nFinalIn, bool isNonDiffIn, bool isResolvedIn, bool isDiffractiveAin, bool isDiffractiveBin, bool isDiffractiveCin)
{
    setType(nameIn, codeIn, nFinalIn, isNonDiffIn, isResolvedIn, isDiffractiveAin, isDiffractiveBin, isDiffractiveCin);
}


void Pythia8::Info::setType__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int codeIn, int nFinalIn, bool isNonDiffIn, bool isResolvedIn, bool isDiffractiveAin, bool isDiffractiveBin)
{
    setType(nameIn, codeIn, nFinalIn, isNonDiffIn, isResolvedIn, isDiffractiveAin, isDiffractiveBin);
}


void Pythia8::Info::setType__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int codeIn, int nFinalIn, bool isNonDiffIn, bool isResolvedIn, bool isDiffractiveAin)
{
    setType(nameIn, codeIn, nFinalIn, isNonDiffIn, isResolvedIn, isDiffractiveAin);
}


void Pythia8::Info::setType__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int codeIn, int nFinalIn, bool isNonDiffIn, bool isResolvedIn)
{
    setType(nameIn, codeIn, nFinalIn, isNonDiffIn, isResolvedIn);
}


void Pythia8::Info::setType__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int codeIn, int nFinalIn, bool isNonDiffIn)
{
    setType(nameIn, codeIn, nFinalIn, isNonDiffIn);
}


void Pythia8::Info::setType__BOSS(std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int codeIn, int nFinalIn)
{
    setType(nameIn, codeIn, nFinalIn);
}


void Pythia8::Info::setTypeMPI__BOSS(int codeMPIIn, double pTMPIIn, int iAMPIIn, int iBMPIIn)
{
    setTypeMPI(codeMPIIn, pTMPIIn, iAMPIIn, iBMPIIn);
}


void Pythia8::Info::setTypeMPI__BOSS(int codeMPIIn, double pTMPIIn, int iAMPIIn)
{
    setTypeMPI(codeMPIIn, pTMPIIn, iAMPIIn);
}


void Pythia8::Info::setTypeMPI__BOSS(int codeMPIIn, double pTMPIIn)
{
    setTypeMPI(codeMPIIn, pTMPIIn);
}


void Pythia8::Info::setImpact__BOSS(double bMPIIn, double enhanceMPIIn)
{
    setImpact(bMPIIn, enhanceMPIIn);
}




#include "backend_types/Pythia_8_186/identification.hpp"

Pythia8::Abstract_Info* Pythia8::Info::pointerCopy__BOSS()
{
    Pythia8::Abstract_Info* new_ptr = new Pythia8::Info(*this);
    new_ptr->can_delete_wrapper(true);
    return new_ptr;
}

void Pythia8::Info::pointerAssign__BOSS(Pythia8::Abstract_Info* in)
{
    CAT_3(BACKENDNAME,_,SAFE_VERSION)::Pythia8::Info* wptr_temp = wrapper__BOSS();
    *this = *dynamic_cast<Info*>(in);
    wrapper__BOSS(wptr_temp);
}

#include "gambit/Backends/backend_undefs.hpp"
