// Forward declarations:
namespace Pythia8
{
    class Abstract__Particle;
    class Particle;
    class Abstract__ParticleData;
    class ParticleData;
    class Abstract__Vec4;
    class Vec4;
    class Abstract__RotBstMatrix;
    class RotBstMatrix;
    class Abstract__Event;
    class Event;
    class Abstract__Pythia;
    class Pythia;
}


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wreturn-type"
namespace Pythia8
{
    class Abstract__Event
    {
        private:
        public:
            // UNKNOWN: OperatorMethod

            virtual void init_GAMBIT(std::string headerIn, Pythia8::Abstract__ParticleData* particleDataPtrIn, int startColTagIn) {};
            void init(std::string headerIn, Pythia8::Abstract__ParticleData* particleDataPtrIn, int startColTagIn)
            {
                init_GAMBIT( headerIn,  particleDataPtrIn,  startColTagIn);
            }

            virtual void clear() {};

            virtual void reset() {};
            // UNKNOWN: OperatorMethod
            // UNKNOWN: OperatorMethod

            virtual Pythia8::Particle*& front_GAMBIT() {};
            Pythia8::Abstract__Particle*& front()
            {
                return reinterpret_cast<Pythia8::Abstract__Particle*&>(front_GAMBIT());
            }

            virtual Pythia8::Particle*& at_GAMBIT(int i) {};
            Pythia8::Abstract__Particle*& at(int i)
            {
                return reinterpret_cast<Pythia8::Abstract__Particle*&>(at_GAMBIT( i));
            }

            virtual Pythia8::Particle*& back_GAMBIT() {};
            Pythia8::Abstract__Particle*& back()
            {
                return reinterpret_cast<Pythia8::Abstract__Particle*&>(back_GAMBIT());
            }

            virtual int size() {};

            virtual int* append_GAMBIT(Pythia8::Abstract__Particle entryIn) {};
            int* append(Pythia8::Abstract__Particle entryIn)
            {
                return append_GAMBIT( entryIn);
            }

            virtual int append(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn, double polIn) {};

            virtual int* append_GAMBIT(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Pythia8::Abstract__Vec4 p, double m, double scaleIn, double polIn) {};
            int* append(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Pythia8::Abstract__Vec4 p, double m, double scaleIn, double polIn)
            {
                return append_GAMBIT( id,  status,  mother1,  mother2,  daughter1,  daughter2,  col,  acol,  p,  m,  scaleIn,  polIn);
            }

            virtual int append(int id, int status, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn, double polIn) {};

            virtual int* append_GAMBIT(int id, int status, int col, int acol, Pythia8::Abstract__Vec4 p, double m, double scaleIn, double polIn) {};
            int* append(int id, int status, int col, int acol, Pythia8::Abstract__Vec4 p, double m, double scaleIn, double polIn)
            {
                return append_GAMBIT( id,  status,  col,  acol,  p,  m,  scaleIn,  polIn);
            }

            virtual void setEvtPtr(int iSet) {};

            virtual int copy(int iCopy, int newStatus) {};

            virtual void list() {};

            virtual void list(std::ostream& os) {};

            virtual void list(bool showScaleAndVertex, bool showMothersAndDaughters) {};

            virtual void list(bool showScaleAndVertex, bool showMothersAndDaughters, std::ostream& os) {};

            virtual void popBack(int nRemove) {};

            virtual void remove(int iFirst, int iLast) {};

            virtual bool undoDecay(int i) {};

            virtual void restorePtrs() {};

            virtual void saveSize() {};

            virtual void restoreSize() {};

            virtual int savedSizeValue() {};

            virtual void initColTag(int colTag) {};

            virtual int lastColTag() {};

            virtual int nextColTag() {};

            virtual void scale(double scaleIn) {};

            virtual double scale() {};

            virtual void scaleSecond(double scaleSecondIn) {};

            virtual double scaleSecond() {};

            virtual int statusHepMC(int i) {};

            virtual int iTopCopy(int i) {};

            virtual int iBotCopy(int i) {};

            virtual int iTopCopyId(int i) {};

            virtual int iBotCopyId(int i) {};

            virtual bool isAncestor(int i, int iAncestor) {};

            virtual void rot(double theta, double phi) {};

            virtual void bst(double betaX, double betaY, double betaZ) {};

            virtual void bst(double betaX, double betaY, double betaZ, double gamma) {};

            virtual void bst_GAMBIT(const Pythia8::Abstract__Vec4& vec) {};
            void bst(const Pythia8::Abstract__Vec4& vec)
            {
                bst_GAMBIT( vec);
            }

            virtual void rotbst_GAMBIT(const Pythia8::Abstract__RotBstMatrix& M) {};
            void rotbst(const Pythia8::Abstract__RotBstMatrix& M)
            {
                rotbst_GAMBIT( M);
            }

            virtual void clearJunctions() {};

            virtual void appendJunction(int kind, int col0, int col1, int col2) {};

            virtual int sizeJunction() {};

            virtual bool remainsJunction(int i) {};

            virtual void remainsJunction(int i, bool remainsIn) {};

            virtual int kindJunction(int i) {};

            virtual int colJunction(int i, int j) {};

            virtual void colJunction(int i, int j, int colIn) {};

            virtual int endColJunction(int i, int j) {};

            virtual void endColJunction(int i, int j, int endColIn) {};

            virtual int statusJunction(int i, int j) {};

            virtual void statusJunction(int i, int j, int statusIn) {};

            virtual void eraseJunction(int i) {};

            virtual void saveJunctionSize() {};

            virtual void restoreJunctionSize() {};

            virtual void listJunctions(std::ostream& os) {};
            // UNKNOWN: OperatorMethod

        public:
            Event* downcast()
            {
                return reinterpret_cast<Event*>(this);
            }
    };
}
#pragma GCC diagnostic pop

