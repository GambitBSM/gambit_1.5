#include <string>
#include "Pythia8/ParticleData.h"
#include "backend_types/Pythia_8_209/wrapper_ParticleDataEntry.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_0(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn, double tau0In)
    {
        return new ParticleDataEntry(idIn, nameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn, mMaxIn, tau0In);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_1(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn)
    {
        return new ParticleDataEntry(idIn, nameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn, mMaxIn);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_2(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn)
    {
        return new ParticleDataEntry(idIn, nameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_3(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn)
    {
        return new ParticleDataEntry(idIn, nameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_4(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In)
    {
        return new ParticleDataEntry(idIn, nameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_5(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn)
    {
        return new ParticleDataEntry(idIn, nameIn, spinTypeIn, chargeTypeIn, colTypeIn);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_6(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int spinTypeIn, int chargeTypeIn)
    {
        return new ParticleDataEntry(idIn, nameIn, spinTypeIn, chargeTypeIn);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_7(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, int spinTypeIn)
    {
        return new ParticleDataEntry(idIn, nameIn, spinTypeIn);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_8(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn)
    {
        return new ParticleDataEntry(idIn, nameIn);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_9(int idIn)
    {
        return new ParticleDataEntry(idIn);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_10()
    {
        return new ParticleDataEntry();
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_11(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn, double tau0In)
    {
        return new ParticleDataEntry(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn, mMaxIn, tau0In);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_12(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn, double mMaxIn)
    {
        return new ParticleDataEntry(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn, mMaxIn);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_13(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn, double mMinIn)
    {
        return new ParticleDataEntry(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn, mMinIn);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_14(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In, double mWidthIn)
    {
        return new ParticleDataEntry(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In, mWidthIn);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_15(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn, double m0In)
    {
        return new ParticleDataEntry(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn, m0In);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_16(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn, int colTypeIn)
    {
        return new ParticleDataEntry(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn, colTypeIn);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_17(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn, int chargeTypeIn)
    {
        return new ParticleDataEntry(idIn, nameIn, antiNameIn, spinTypeIn, chargeTypeIn);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_18(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn, int spinTypeIn)
    {
        return new ParticleDataEntry(idIn, nameIn, antiNameIn, spinTypeIn);
    }
    
    Abstract_ParticleDataEntry* Factory_ParticleDataEntry_19(int idIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > nameIn, std::basic_string<char,std::char_traits<char>,std::allocator<char> > antiNameIn)
    {
        return new ParticleDataEntry(idIn, nameIn, antiNameIn);
    }
    
}

