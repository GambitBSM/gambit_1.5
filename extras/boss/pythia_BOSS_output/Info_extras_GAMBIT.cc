#include <map>
#include <string>
#include <vector>
#include <ostream>
#include "abstracts_typedefs.hpp"
#include "wrappers_typedefs.hpp"
#include "Pythia8/Info.h"

void Pythia8::Info::list_GAMBIT() const
{
    list();
}


bool Pythia8::Info::hasSub_GAMBIT() const
{
    return hasSub();
}


std::string Pythia8::Info::nameSub_GAMBIT() const
{
    return nameSub();
}


int Pythia8::Info::codeSub_GAMBIT() const
{
    return codeSub();
}


int Pythia8::Info::nFinalSub_GAMBIT() const
{
    return nFinalSub();
}


int Pythia8::Info::id1_GAMBIT() const
{
    return id1();
}


int Pythia8::Info::id2_GAMBIT() const
{
    return id2();
}


double Pythia8::Info::x1_GAMBIT() const
{
    return x1();
}


double Pythia8::Info::x2_GAMBIT() const
{
    return x2();
}


double Pythia8::Info::y_GAMBIT() const
{
    return y();
}


double Pythia8::Info::tau_GAMBIT() const
{
    return tau();
}


int Pythia8::Info::id1pdf_GAMBIT() const
{
    return id1pdf();
}


int Pythia8::Info::id2pdf_GAMBIT() const
{
    return id2pdf();
}


double Pythia8::Info::x1pdf_GAMBIT() const
{
    return x1pdf();
}


double Pythia8::Info::x2pdf_GAMBIT() const
{
    return x2pdf();
}


double Pythia8::Info::pdf1_GAMBIT() const
{
    return pdf1();
}


double Pythia8::Info::pdf2_GAMBIT() const
{
    return pdf2();
}


double Pythia8::Info::QFac_GAMBIT() const
{
    return QFac();
}


double Pythia8::Info::Q2Fac_GAMBIT() const
{
    return Q2Fac();
}


double Pythia8::Info::alphaS_GAMBIT() const
{
    return alphaS();
}


double Pythia8::Info::alphaEM_GAMBIT() const
{
    return alphaEM();
}


double Pythia8::Info::QRen_GAMBIT() const
{
    return QRen();
}


double Pythia8::Info::Q2Ren_GAMBIT() const
{
    return Q2Ren();
}


double Pythia8::Info::scalup_GAMBIT() const
{
    return scalup();
}


double Pythia8::Info::mHat_GAMBIT() const
{
    return mHat();
}


double Pythia8::Info::sHat_GAMBIT() const
{
    return sHat();
}


double Pythia8::Info::tHat_GAMBIT() const
{
    return tHat();
}


double Pythia8::Info::uHat_GAMBIT() const
{
    return uHat();
}


double Pythia8::Info::pTHat_GAMBIT() const
{
    return pTHat();
}


double Pythia8::Info::pT2Hat_GAMBIT() const
{
    return pT2Hat();
}


double Pythia8::Info::m3Hat_GAMBIT() const
{
    return m3Hat();
}


double Pythia8::Info::m4Hat_GAMBIT() const
{
    return m4Hat();
}


double Pythia8::Info::thetaHat_GAMBIT() const
{
    return thetaHat();
}


double Pythia8::Info::phiHat_GAMBIT() const
{
    return phiHat();
}


std::string Pythia8::Info::nameProc_GAMBIT()
{
    return nameProc();
}


long int Pythia8::Info::nTried_GAMBIT()
{
    return nTried();
}


long int Pythia8::Info::nSelected_GAMBIT()
{
    return nSelected();
}


long int Pythia8::Info::nAccepted_GAMBIT()
{
    return nAccepted();
}


double Pythia8::Info::sigmaGen_GAMBIT()
{
    return sigmaGen();
}


double Pythia8::Info::sigmaErr_GAMBIT()
{
    return sigmaErr();
}


void Pythia8::Info::setCounter_GAMBIT(int i)
{
    setCounter(i);
}


void Pythia8::Info::addCounter_GAMBIT(int i)
{
    addCounter(i);
}


void Pythia8::Info::errorMsg_GAMBIT(std::string messageIn, std::string extraIn, bool showAlways)
{
    errorMsg(messageIn, extraIn, showAlways);
}


void Pythia8::Info::errorMsg_GAMBIT(std::string messageIn, std::string extraIn)
{
    errorMsg(messageIn, extraIn);
}


void Pythia8::Info::errorMsg_GAMBIT(std::string messageIn)
{
    errorMsg(messageIn);
}


void Pythia8::Info::errorStatistics_GAMBIT()
{
    errorStatistics();
}


void Pythia8::Info::setType_GAMBIT(std::string nameIn, int codeIn, int nFinalIn, bool isNonDiffIn, bool isResolvedIn, bool isDiffractiveAin, bool isDiffractiveBin, bool isDiffractiveCin)
{
    setType(nameIn, codeIn, nFinalIn, isNonDiffIn, isResolvedIn, isDiffractiveAin, isDiffractiveBin, isDiffractiveCin);
}


void Pythia8::Info::setType_GAMBIT(std::string nameIn, int codeIn, int nFinalIn, bool isNonDiffIn, bool isResolvedIn, bool isDiffractiveAin, bool isDiffractiveBin)
{
    setType(nameIn, codeIn, nFinalIn, isNonDiffIn, isResolvedIn, isDiffractiveAin, isDiffractiveBin);
}


void Pythia8::Info::setType_GAMBIT(std::string nameIn, int codeIn, int nFinalIn, bool isNonDiffIn, bool isResolvedIn, bool isDiffractiveAin)
{
    setType(nameIn, codeIn, nFinalIn, isNonDiffIn, isResolvedIn, isDiffractiveAin);
}


void Pythia8::Info::setType_GAMBIT(std::string nameIn, int codeIn, int nFinalIn, bool isNonDiffIn, bool isResolvedIn)
{
    setType(nameIn, codeIn, nFinalIn, isNonDiffIn, isResolvedIn);
}


void Pythia8::Info::setType_GAMBIT(std::string nameIn, int codeIn, int nFinalIn, bool isNonDiffIn)
{
    setType(nameIn, codeIn, nFinalIn, isNonDiffIn);
}


void Pythia8::Info::setType_GAMBIT(std::string nameIn, int codeIn, int nFinalIn)
{
    setType(nameIn, codeIn, nFinalIn);
}


void Pythia8::Info::setTypeMPI_GAMBIT(int codeMPIIn, double pTMPIIn, int iAMPIIn, int iBMPIIn)
{
    setTypeMPI(codeMPIIn, pTMPIIn, iAMPIIn, iBMPIIn);
}


void Pythia8::Info::setTypeMPI_GAMBIT(int codeMPIIn, double pTMPIIn, int iAMPIIn)
{
    setTypeMPI(codeMPIIn, pTMPIIn, iAMPIIn);
}


void Pythia8::Info::setTypeMPI_GAMBIT(int codeMPIIn, double pTMPIIn)
{
    setTypeMPI(codeMPIIn, pTMPIIn);
}


void Pythia8::Info::setImpact_GAMBIT(double bMPIIn, double enhanceMPIIn)
{
    setImpact(bMPIIn, enhanceMPIIn);
}



Pythia8::Abstract_Info* Pythia8::Info::pointerCopy_GAMBIT() { return new Pythia8::Info(*this); }
void Pythia8::Info::pointerAssign_GAMBIT(Pythia8::Abstract_Info* in) { *this = *dynamic_cast<Info*>(in); }
