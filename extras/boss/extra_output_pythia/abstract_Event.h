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

            virtual void clear_GAMBIT() {};
            void clear()
            {
                clear_GAMBIT();
            }

            virtual void reset_GAMBIT() {};
            void reset()
            {
                reset_GAMBIT();
            }
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

            virtual int* size_GAMBIT() {};
            int* size()
            {
                return size_GAMBIT();
            }

            virtual int* append_GAMBIT(Pythia8::Abstract__Particle entryIn) {};
            int* append(Pythia8::Abstract__Particle entryIn)
            {
                return append_GAMBIT( entryIn);
            }

            virtual int* append_GAMBIT(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn, double polIn) {};
            int* append(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn, double polIn)
            {
                return append_GAMBIT( id,  status,  mother1,  mother2,  daughter1,  daughter2,  col,  acol,  px,  py,  pz,  e,  m,  scaleIn,  polIn);
            }

            virtual int* append_GAMBIT(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Pythia8::Abstract__Vec4 p, double m, double scaleIn, double polIn) {};
            int* append(int id, int status, int mother1, int mother2, int daughter1, int daughter2, int col, int acol, Pythia8::Abstract__Vec4 p, double m, double scaleIn, double polIn)
            {
                return append_GAMBIT( id,  status,  mother1,  mother2,  daughter1,  daughter2,  col,  acol,  p,  m,  scaleIn,  polIn);
            }

            virtual int* append_GAMBIT(int id, int status, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn, double polIn) {};
            int* append(int id, int status, int col, int acol, double px, double py, double pz, double e, double m, double scaleIn, double polIn)
            {
                return append_GAMBIT( id,  status,  col,  acol,  px,  py,  pz,  e,  m,  scaleIn,  polIn);
            }

            virtual int* append_GAMBIT(int id, int status, int col, int acol, Pythia8::Abstract__Vec4 p, double m, double scaleIn, double polIn) {};
            int* append(int id, int status, int col, int acol, Pythia8::Abstract__Vec4 p, double m, double scaleIn, double polIn)
            {
                return append_GAMBIT( id,  status,  col,  acol,  p,  m,  scaleIn,  polIn);
            }

            virtual void setEvtPtr_GAMBIT(int iSet) {};
            void setEvtPtr(int iSet)
            {
                setEvtPtr_GAMBIT( iSet);
            }

            virtual int* copy_GAMBIT(int iCopy, int newStatus) {};
            int* copy(int iCopy, int newStatus)
            {
                return copy_GAMBIT( iCopy,  newStatus);
            }

            virtual void list_GAMBIT() {};
            void list()
            {
                list_GAMBIT();
            }

            virtual void list_GAMBIT(std::ostream& os) {};
            void list(std::ostream& os)
            {
                list_GAMBIT( os);
            }

            virtual void list_GAMBIT(bool showScaleAndVertex, bool showMothersAndDaughters) {};
            void list(bool showScaleAndVertex, bool showMothersAndDaughters)
            {
                list_GAMBIT( showScaleAndVertex,  showMothersAndDaughters);
            }

            virtual void list_GAMBIT(bool showScaleAndVertex, bool showMothersAndDaughters, std::ostream& os) {};
            void list(bool showScaleAndVertex, bool showMothersAndDaughters, std::ostream& os)
            {
                list_GAMBIT( showScaleAndVertex,  showMothersAndDaughters,  os);
            }

            virtual void popBack_GAMBIT(int nRemove) {};
            void popBack(int nRemove)
            {
                popBack_GAMBIT( nRemove);
            }

            virtual void remove_GAMBIT(int iFirst, int iLast) {};
            void remove(int iFirst, int iLast)
            {
                remove_GAMBIT( iFirst,  iLast);
            }

            virtual bool* undoDecay_GAMBIT(int i) {};
            bool* undoDecay(int i)
            {
                return undoDecay_GAMBIT( i);
            }

            virtual void restorePtrs_GAMBIT() {};
            void restorePtrs()
            {
                restorePtrs_GAMBIT();
            }

            virtual void saveSize_GAMBIT() {};
            void saveSize()
            {
                saveSize_GAMBIT();
            }

            virtual void restoreSize_GAMBIT() {};
            void restoreSize()
            {
                restoreSize_GAMBIT();
            }

            virtual int* savedSizeValue_GAMBIT() {};
            int* savedSizeValue()
            {
                return savedSizeValue_GAMBIT();
            }

            virtual void initColTag_GAMBIT(int colTag) {};
            void initColTag(int colTag)
            {
                initColTag_GAMBIT( colTag);
            }

            virtual int* lastColTag_GAMBIT() {};
            int* lastColTag()
            {
                return lastColTag_GAMBIT();
            }

            virtual int* nextColTag_GAMBIT() {};
            int* nextColTag()
            {
                return nextColTag_GAMBIT();
            }

            virtual void scale_GAMBIT(double scaleIn) {};
            void scale(double scaleIn)
            {
                scale_GAMBIT( scaleIn);
            }

            virtual double* scale_GAMBIT() {};
            double* scale()
            {
                return scale_GAMBIT();
            }

            virtual void scaleSecond_GAMBIT(double scaleSecondIn) {};
            void scaleSecond(double scaleSecondIn)
            {
                scaleSecond_GAMBIT( scaleSecondIn);
            }

            virtual double* scaleSecond_GAMBIT() {};
            double* scaleSecond()
            {
                return scaleSecond_GAMBIT();
            }

            virtual int* statusHepMC_GAMBIT(int i) {};
            int* statusHepMC(int i)
            {
                return statusHepMC_GAMBIT( i);
            }

            virtual int* iTopCopy_GAMBIT(int i) {};
            int* iTopCopy(int i)
            {
                return iTopCopy_GAMBIT( i);
            }

            virtual int* iBotCopy_GAMBIT(int i) {};
            int* iBotCopy(int i)
            {
                return iBotCopy_GAMBIT( i);
            }

            virtual int* iTopCopyId_GAMBIT(int i) {};
            int* iTopCopyId(int i)
            {
                return iTopCopyId_GAMBIT( i);
            }

            virtual int* iBotCopyId_GAMBIT(int i) {};
            int* iBotCopyId(int i)
            {
                return iBotCopyId_GAMBIT( i);
            }

            virtual bool* isAncestor_GAMBIT(int i, int iAncestor) {};
            bool* isAncestor(int i, int iAncestor)
            {
                return isAncestor_GAMBIT( i,  iAncestor);
            }

            virtual void rot_GAMBIT(double theta, double phi) {};
            void rot(double theta, double phi)
            {
                rot_GAMBIT( theta,  phi);
            }

            virtual void bst_GAMBIT(double betaX, double betaY, double betaZ) {};
            void bst(double betaX, double betaY, double betaZ)
            {
                bst_GAMBIT( betaX,  betaY,  betaZ);
            }

            virtual void bst_GAMBIT(double betaX, double betaY, double betaZ, double gamma) {};
            void bst(double betaX, double betaY, double betaZ, double gamma)
            {
                bst_GAMBIT( betaX,  betaY,  betaZ,  gamma);
            }

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

            virtual void clearJunctions_GAMBIT() {};
            void clearJunctions()
            {
                clearJunctions_GAMBIT();
            }

            virtual void appendJunction_GAMBIT(int kind, int col0, int col1, int col2) {};
            void appendJunction(int kind, int col0, int col1, int col2)
            {
                appendJunction_GAMBIT( kind,  col0,  col1,  col2);
            }

            virtual int* sizeJunction_GAMBIT() {};
            int* sizeJunction()
            {
                return sizeJunction_GAMBIT();
            }

            virtual bool* remainsJunction_GAMBIT(int i) {};
            bool* remainsJunction(int i)
            {
                return remainsJunction_GAMBIT( i);
            }

            virtual void remainsJunction_GAMBIT(int i, bool remainsIn) {};
            void remainsJunction(int i, bool remainsIn)
            {
                remainsJunction_GAMBIT( i,  remainsIn);
            }

            virtual int* kindJunction_GAMBIT(int i) {};
            int* kindJunction(int i)
            {
                return kindJunction_GAMBIT( i);
            }

            virtual int* colJunction_GAMBIT(int i, int j) {};
            int* colJunction(int i, int j)
            {
                return colJunction_GAMBIT( i,  j);
            }

            virtual void colJunction_GAMBIT(int i, int j, int colIn) {};
            void colJunction(int i, int j, int colIn)
            {
                colJunction_GAMBIT( i,  j,  colIn);
            }

            virtual int* endColJunction_GAMBIT(int i, int j) {};
            int* endColJunction(int i, int j)
            {
                return endColJunction_GAMBIT( i,  j);
            }

            virtual void endColJunction_GAMBIT(int i, int j, int endColIn) {};
            void endColJunction(int i, int j, int endColIn)
            {
                endColJunction_GAMBIT( i,  j,  endColIn);
            }

            virtual int* statusJunction_GAMBIT(int i, int j) {};
            int* statusJunction(int i, int j)
            {
                return statusJunction_GAMBIT( i,  j);
            }

            virtual void statusJunction_GAMBIT(int i, int j, int statusIn) {};
            void statusJunction(int i, int j, int statusIn)
            {
                statusJunction_GAMBIT( i,  j,  statusIn);
            }

            virtual void eraseJunction_GAMBIT(int i) {};
            void eraseJunction(int i)
            {
                eraseJunction_GAMBIT( i);
            }

            virtual void saveJunctionSize_GAMBIT() {};
            void saveJunctionSize()
            {
                saveJunctionSize_GAMBIT();
            }

            virtual void restoreJunctionSize_GAMBIT() {};
            void restoreJunctionSize()
            {
                restoreJunctionSize_GAMBIT();
            }

            virtual void listJunctions_GAMBIT(std::ostream& os) {};
            void listJunctions(std::ostream& os)
            {
                listJunctions_GAMBIT( os);
            }
            // UNKNOWN: OperatorMethod

        public:
            Event* downcast()
            {
                return reinterpret_cast<Event*>(this);
            }
    };
}

