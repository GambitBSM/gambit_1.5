#include "Pythia8/Event.h"
#include "backend_types/Pythia_8_186/wrapper_Particle.h"
#include "Pythia8/Basics.h"
#include "backend_types/Pythia_8_186/wrapper_Vec4.h"
#include "gambit/Backends/abstracttypedefs.h"
#include "gambit/Backends/wrappertypedefs.h"

namespace Pythia8
{
    Abstract_Particle* Factory_Particle_0()
    {
        return new Particle();
    }
    
    Abstract_Particle* Factory_Particle_1(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn, double pxIn, double pyIn, double pzIn, double eIn, double mIn, double scaleIn, double polIn)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In, colIn, acolIn, pxIn, pyIn, pzIn, eIn, mIn, scaleIn, polIn);
    }
    
    Abstract_Particle* Factory_Particle_2(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn, double pxIn, double pyIn, double pzIn, double eIn, double mIn, double scaleIn)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In, colIn, acolIn, pxIn, pyIn, pzIn, eIn, mIn, scaleIn);
    }
    
    Abstract_Particle* Factory_Particle_3(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn, double pxIn, double pyIn, double pzIn, double eIn, double mIn)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In, colIn, acolIn, pxIn, pyIn, pzIn, eIn, mIn);
    }
    
    Abstract_Particle* Factory_Particle_4(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn, double pxIn, double pyIn, double pzIn, double eIn)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In, colIn, acolIn, pxIn, pyIn, pzIn, eIn);
    }
    
    Abstract_Particle* Factory_Particle_5(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn, double pxIn, double pyIn, double pzIn)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In, colIn, acolIn, pxIn, pyIn, pzIn);
    }
    
    Abstract_Particle* Factory_Particle_6(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn, double pxIn, double pyIn)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In, colIn, acolIn, pxIn, pyIn);
    }
    
    Abstract_Particle* Factory_Particle_7(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn, double pxIn)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In, colIn, acolIn, pxIn);
    }
    
    Abstract_Particle* Factory_Particle_8(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In, colIn, acolIn);
    }
    
    Abstract_Particle* Factory_Particle_9(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In, colIn);
    }
    
    Abstract_Particle* Factory_Particle_10(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In);
    }
    
    Abstract_Particle* Factory_Particle_11(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In, daughter1In);
    }
    
    Abstract_Particle* Factory_Particle_12(int idIn, int statusIn, int mother1In, int mother2In)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In);
    }
    
    Abstract_Particle* Factory_Particle_13(int idIn, int statusIn, int mother1In)
    {
        return new Particle(idIn, statusIn, mother1In);
    }
    
    Abstract_Particle* Factory_Particle_14(int idIn, int statusIn)
    {
        return new Particle(idIn, statusIn);
    }
    
    Abstract_Particle* Factory_Particle_15(int idIn)
    {
        return new Particle(idIn);
    }
    
    Abstract_Particle* Factory_Particle_16(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn, Pythia8::Vec4__BOSS& pIn, double mIn, double scaleIn, double polIn)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In, colIn, acolIn, dynamic_cast< Pythia8::Vec4& >(*pIn.BEptr), mIn, scaleIn, polIn);
    }
    
    Abstract_Particle* Factory_Particle_17(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn, Pythia8::Vec4__BOSS& pIn, double mIn, double scaleIn)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In, colIn, acolIn, dynamic_cast< Pythia8::Vec4& >(*pIn.BEptr), mIn, scaleIn);
    }
    
    Abstract_Particle* Factory_Particle_18(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn, Pythia8::Vec4__BOSS& pIn, double mIn)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In, colIn, acolIn, dynamic_cast< Pythia8::Vec4& >(*pIn.BEptr), mIn);
    }
    
    Abstract_Particle* Factory_Particle_19(int idIn, int statusIn, int mother1In, int mother2In, int daughter1In, int daughter2In, int colIn, int acolIn, Pythia8::Vec4__BOSS& pIn)
    {
        return new Particle(idIn, statusIn, mother1In, mother2In, daughter1In, daughter2In, colIn, acolIn, dynamic_cast< Pythia8::Vec4& >(*pIn.BEptr));
    }
    
}

