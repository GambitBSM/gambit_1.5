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
    class Abstract__RotBstMatrix
    {
        private:
        public:
            // UNKNOWN: OperatorMethod

            virtual void rot_GAMBIT(double arg_1, double arg_2) {};
            void rot(double arg_1, double arg_2)
            {
                rot_GAMBIT( arg_1,  arg_2);
            }

            virtual void rot_GAMBIT(const Pythia8::Abstract__Vec4& p) {};
            void rot(const Pythia8::Abstract__Vec4& p)
            {
                rot_GAMBIT( p);
            }

            virtual void bst_GAMBIT(double arg_1, double arg_2, double arg_3) {};
            void bst(double arg_1, double arg_2, double arg_3)
            {
                bst_GAMBIT( arg_1,  arg_2,  arg_3);
            }

            virtual void bst_GAMBIT(const Pythia8::Abstract__Vec4& arg_1) {};
            void bst(const Pythia8::Abstract__Vec4& arg_1)
            {
                bst_GAMBIT( arg_1);
            }

            virtual void bstback_GAMBIT(const Pythia8::Abstract__Vec4& arg_1) {};
            void bstback(const Pythia8::Abstract__Vec4& arg_1)
            {
                bstback_GAMBIT( arg_1);
            }

            virtual void bst_GAMBIT(const Pythia8::Abstract__Vec4& arg_1, const Pythia8::Abstract__Vec4& arg_2) {};
            void bst(const Pythia8::Abstract__Vec4& arg_1, const Pythia8::Abstract__Vec4& arg_2)
            {
                bst_GAMBIT( arg_1,  arg_2);
            }

            virtual void toCMframe_GAMBIT(const Pythia8::Abstract__Vec4& arg_1, const Pythia8::Abstract__Vec4& arg_2) {};
            void toCMframe(const Pythia8::Abstract__Vec4& arg_1, const Pythia8::Abstract__Vec4& arg_2)
            {
                toCMframe_GAMBIT( arg_1,  arg_2);
            }

            virtual void fromCMframe_GAMBIT(const Pythia8::Abstract__Vec4& arg_1, const Pythia8::Abstract__Vec4& arg_2) {};
            void fromCMframe(const Pythia8::Abstract__Vec4& arg_1, const Pythia8::Abstract__Vec4& arg_2)
            {
                fromCMframe_GAMBIT( arg_1,  arg_2);
            }

            virtual void rotbst_GAMBIT(const Pythia8::Abstract__RotBstMatrix& arg_1) {};
            void rotbst(const Pythia8::Abstract__RotBstMatrix& arg_1)
            {
                rotbst_GAMBIT( arg_1);
            }

            virtual void invert_GAMBIT() {};
            void invert()
            {
                invert_GAMBIT();
            }

            virtual void reset_GAMBIT() {};
            void reset()
            {
                reset_GAMBIT();
            }

            virtual double* deviation_GAMBIT() {};
            double* deviation()
            {
                return deviation_GAMBIT();
            }

        public:
            RotBstMatrix* downcast()
            {
                return reinterpret_cast<RotBstMatrix*>(this);
            }
    };
}
#pragma GCC diagnostic pop

