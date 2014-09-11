#include <iostream> 
#include <algorithm>
#include <vector>
#include <math.h>
#include <map>
#include <fstream>

using std::vector;
using std::cout;
using std::endl;
using std::map;
using std::ofstream;
using std::ios;

// Forward declarations
struct vec3;
struct vec4;
struct mat4;
vec4 Ep4vec(const vec3 p, double m);

// ****************************************************************
// Some minimal vector and matrix classes needed in the calculation
// ****************************************************************

// 3-vector
struct vec3
{
    double vals[3];
    // Constructors
    vec3(){}
    vec3(double v0, double v1, double v2){vals[0]=v0;vals[1]=v1;vals[2]=v2;}
    vec3(double v0){vals[0]=v0;vals[1]=v0;vals[2]=v0;}    
    // Get the length of the vector
    double length() const{return sqrt(vals[0]*vals[0]+vals[1]*vals[1]+vals[2]*vals[2]);}
    // Normalize vector to unit length
    void normalize()
    {
        double l = length();
        for(int i=0;i<3;i++)
        {
            vals[i] /= l;
        }        
    }
    // Normalize vector to length len
    void normalize(double len)
    {
        double l = length();
        for(int i=0;i<3;i++)
        {
            vals[i] *= len/l;
        }        
    }    
    // Operators
    double& operator[](unsigned int i){ return vals[i];}
	const double operator[](unsigned int i)const{ return vals[i];}      
    vec3 operator-()const{return vec3(-vals[0],-vals[1],-vals[2]);	} 
	// Non-member operators (not strictly necessary to list them here since all members are public)
    friend vec3 operator* (double x, const vec3 &y);
    friend vec3 operator* (const vec3 &y, double x);
    friend vec3 operator/ (const vec3 &y, double x);	
    friend std::ostream& operator<<(std::ostream& os, const vec3& v);
};
vec3 operator* (double x, const vec3 &y)
{
    vec3 retVec = y;
    for(int i=0;i<3;i++)
    {
        retVec[i] *= x;
    }
    return retVec;
}
vec3 operator* (const vec3 &y, double x)
{
    vec3 retVec = y;
    for(int i=0;i<3;i++)
    {
        retVec[i] *= x;
    }
    return retVec;
}
vec3 operator/ (const vec3 &y, double x)
{
    vec3 retVec = y;
    for(int i=0;i<3;i++)
    {
        retVec[i] /= x;
    }
    return retVec;
}
std::ostream& operator<<(std::ostream& os, const vec3& v)
{
    os << v[0] << ", " << v[1]  << ", " << v[2];
    return os;
}

// 4-vector
struct vec4
{
    double vals[4];     
    // Constructors
    vec4(){}
    vec4(double v0, double v1, double v2, double v3){vals[0]=v0;vals[1]=v1;vals[2]=v2;vals[3]=v3;} 
    vec4(double v0){vals[0]=v0;vals[1]=v0;vals[2]=v0;vals[3]=v0;}   
    vec4(double v0, vec3 v){vals[0]=v0;vals[1]=v[0];vals[2]=v[1];vals[3]=v[2];}   
    // Returns vec3 containting elements 1,2,3
    vec3 xyz() const{return vec3(vals[1],vals[2],vals[3]);}   
    // Operators 
    double& operator[](unsigned int i){ return vals[i];}     
	const double operator[](unsigned int i)const{ return vals[i];}    
	vec4 operator-(){return vec4(-vals[0],-vals[1],-vals[2],-vals[3]);}
	// Non-member operators (not strictly necessary to list them here since all members are public)	
	friend vec4 operator* (double x, const vec4 &y);
	friend vec4 operator* (const vec4 &y, double x);
    friend vec4 operator+ (const vec4 &x, const vec4 &y);
    friend vec4 operator- (const vec4 &x, const vec4 &y);  
	friend std::ostream& operator<<(std::ostream& os, const vec4& v);
	friend vec4 operator* (const mat4 &m, const vec4 &v);
};
vec4 operator* (double x, const vec4 &y)
{
    vec4 retVec = y;
    for(int i=0;i<4;i++)
    {
        retVec[i] *= x;
    }
    return retVec;
}
vec4 operator* (const vec4 &y, double x)
{
    vec4 retVec = y;
    for(int i=0;i<4;i++)
    {
        retVec[i] *= x;
    }
    return retVec;
}
vec4 operator+ (const vec4 &x, const vec4 &y)
{
    vec4 retVec = x;
    for(int i=0;i<4;i++)
    {
        retVec[i] += y[i];
    }
    return retVec;
}
vec4 operator- (const vec4 &x, const vec4 &y)
{
    vec4 retVec = x;
    for(int i=0;i<4;i++)
    {
        retVec[i] -= y[i];
    }
    return retVec;
}
std::ostream& operator<<(std::ostream& os, const vec4& v)
{
    os << v[0] << ", " << v[1]  << ", " << v[2] << ", " << v[3];
    return os;
}

// 4x4 matrix
struct mat4
{
    double vals[4][4];
    mat4(){};
    mat4(   double v00, double v01, double v02, double v03,
            double v10, double v11, double v12, double v13,
            double v20, double v21, double v22, double v23,
            double v30, double v31, double v32, double v33)
    {
        vals[0][0] = v00;
        vals[0][1] = v01;
        vals[0][2] = v02;
        vals[0][3] = v03; 
        vals[1][0] = v10;
        vals[1][1] = v11;
        vals[1][2] = v12;
        vals[1][3] = v13;
        vals[2][0] = v20;
        vals[2][1] = v21;
        vals[2][2] = v22;
        vals[2][3] = v23;
        vals[3][0] = v30;
        vals[3][1] = v31;
        vals[3][2] = v32;
        vals[3][3] = v33;
    }
    mat4(double v)
    {
        vals[0][0] = v;
        vals[0][1] = v;
        vals[0][2] = v;
        vals[0][3] = v; 
        vals[1][0] = v;
        vals[1][1] = v;
        vals[1][2] = v;
        vals[1][3] = v;
        vals[2][0] = v;
        vals[2][1] = v;
        vals[2][2] = v;
        vals[2][3] = v;
        vals[3][0] = v;
        vals[3][1] = v;
        vals[3][2] = v;
        vals[3][3] = v;
    } 
    // Identity matrix   
    static mat4 identity()
    {
        return mat4(1,0,0,0,
                    0,1,0,0,
                    0,0,1,0,
                    0,0,0,1);
    }    
	// Non-member operators (not strictly necessary to list them here since all members are public)    
    friend vec4 operator* (const mat4 &m, const vec4 &v);
    friend mat4 operator* (const mat4 &m1, const mat4 &m2);
    friend std::ostream& operator<<(std::ostream& os, const mat4& m);
};
vec4 operator* (const mat4 &m, const vec4 &v)
{
    vec4 out(0);
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            out[i] += m.vals[i][j]*v[j];
        }
    }
    return out;
}
mat4 operator* (const mat4 &m1, const mat4 &m2)
{
    mat4 out(0);
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            for(int k=0; k<4; k++)
            {
                out.vals[i][j] += m1.vals[i][k]*m2.vals[k][j];
            }
        }
    }
    return out;
}
std::ostream& operator<<(std::ostream& os, const mat4& m)
{
    for(int i=0;i<4;i++)
    {
        for(int j=0;j<4;j++)
        {
            os << m.vals[i][j] << ", ";
        }
        if(i!=3)
            os << endl;
    }
    return os;
}
double dot(const vec3 &a, const vec3 &b)
{
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
double dot(const vec4 &a, const vec4 &b)
{
    return a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3];
}

// Generate a random number between -1 and 1
double rand_m1_1()
{
    return -1.0 + 2.0 * static_cast<double>(rand()) / RAND_MAX;
}

// Generate a random number between 0 and 1
double rand_0_1()
{
    return static_cast<double>(rand()) / RAND_MAX;
}

// Generate a 3-vector to a random point on the unit sphere
vec3 randOnSphere()
{
    double r1,r2;
    do
    {
        r1 = rand_m1_1();
        r2 = rand_m1_1();
    }
    while(r1*r1+r2*r2 >=1.0);
    vec3 v;
    v[0] = 2.0*r1*sqrt(1-r1*r1-r2*r2);
    v[1] = 2.0*r2*sqrt(1-r1*r1-r2*r2);
    v[2] = 1.0-2.0*(r1*r1+r2*r2);
    return v;
}

// Calculate Lorentz boost matrix corresponding to beta_xyz
void lorentzMatrix(const vec3 &beta_xyz, mat4 &mat)
{
    double b = beta_xyz.length();
    double bm2 = b==0 ? 0 : 1.0/(b*b);
    double bx = beta_xyz[0];
    double by = beta_xyz[1];
    double bz = beta_xyz[2];                 
    double g = 1.0/sqrt(1-b*b);
    mat =  mat4(    g,      -g*bx,              -g*by,              -g*bz,
                    -g*bx,  1+(g-1)*bx*bx*bm2,  (g-1)*bx*by*bm2,    (g-1)*bx*bz*bm2,
                    -g*by,  (g-1)*by*bx*bm2,    1+(g-1)*by*by*bm2,  (g-1)*by*bz*bm2,
                    -g*bz,  (g-1)*bz*bx*bm2,    (g-1)*bz*by*bm2,    1+(g-1)*bz*bz*bm2);  
}
// Boost inVec according to beta_xyz
vec4 lorentzBoost(const vec4 &inVec, const vec3 &beta_xyz)
{
    mat4 lorentz;
    lorentzMatrix(beta_xyz, lorentz);
    return lorentz*inVec;
}

// Boost inVec to the frame where a particle at rest (in this frame) would have 4-momentum p_parent
vec4 p_parentFrame(const vec4 &inVec, const vec4 &p_parent)
{
    vec3 beta_xyz = -p_parent.xyz()/p_parent[0];
    return lorentzBoost(inVec, beta_xyz);
}

void boostMatrixParentFrame(mat4 &mat, vec4 &p_parent)
{
    vec3 beta_xyz = -p_parent.xyz()/p_parent[0];     
    lorentzMatrix(beta_xyz, mat);
}

// Decay channel. Contains the partial width of the channel and the PIDs of the final state particles
struct channel
{
    channel(int ipID0, int ipID1, double ipartialWidth)
    {
        pID0 = ipID0;
        pID1 = ipID1;
        partialWidth = ipartialWidth;
    }
    int pID0;
    int pID1;
    double partialWidth;
    bool same(const channel &in) const{return in.pID0==pID0 && in.pID1 == pID1;}
	int& operator[](unsigned int i)
    {
        if(i<2)
            return (i==0) ? pID0 : pID1;
        else
        {
            cout << "Error: Trying to access non-existing decay final state." << endl;
            exit(1);
            return pID0;
        }            
    }            
};

// Class containing all the allowed decay channels of a single particle, as well as the particle mass.
class DecayTableEntry
{
    public:
        const double m;
        bool chainDecay;        
        DecayTableEntry(int iPID, double im, bool ichainDecay) : m(im)
        {
            pID = iPID;
            chainDecay = ichainDecay;
            totalWidth = 0;
            randInit = false;
            Nchn = 0;
        }
        DecayTableEntry() : m(-1)
        {
            pID = 0;
            chainDecay = false;
            totalWidth = 0;
            randInit = false;
            Nchn = 0;
        }      
        int findChannelIdx(double pick)
        {
            if(!randInit)
                generateRandTable();            
            vector<double>::iterator pos = upper_bound(randLims.begin(),randLims.end(),pick);   
            return pos - randLims.begin();
        }        
        channel randomDecay()
        {
            double pick = rand_0_1();
            int idx = findChannelIdx(pick);
            return table[idx];
        }    
        void generateRandTable()
        {
            randLims.clear();            
            double tmp=0;
            for(int i=0;i<table.size();i++)
            {
                tmp += table[i].partialWidth/totalWidth;
                randLims.push_back(tmp);
            } 
            randInit=true;           
        }
        bool allowedDecay(const channel &in) const
        {
            for(int i=0;i<table.size();i++)
            {
                if(table[i].same(in)) return true;
            }
            return false;
        }
        void addChannel(int pID0, int pID1, double partialWidth)
        {
            if(!chainDecay)
            {
                cout << "Warning: Adding decay channel to particle that is not chain decayed." << endl
                     << "Channel will be ignored unless particle is explicitly set to decay."  << endl;
            }            
            channel chn(pID0, pID1, partialWidth);
            table.push_back(chn);
            totalWidth += partialWidth;
            Nchn++;
        }
        void addChannel(const channel &in)
        {
            if(!chainDecay)
            {
                cout << "Warning: Adding decay channel to particle that is not chain decayed." << endl
                     << "Channel will be ignored unless chainDecay is set to true."  << endl;
            }            
            table.push_back(in);
            totalWidth += in.partialWidth;
            Nchn++;
        }       
        double get_m(){return m;}
        unsigned int get_Nchn(){return Nchn;}
    private:
        vector<channel> table;
        vector<double> randLims;
        double totalWidth; 
        int pID;
        bool randInit;
        unsigned int Nchn;
};

// Table of all particles and their decay channels.
// Uses particle PID as array index
struct DecayTable
{
    map<int,DecayTableEntry> table;
    // Add particle to decay table, specifying particle ID, mass and whether or not it should be decayed in decay chains
    void addEntry(int iPID, double im, bool ichainDecay)
    {
        table.insert ( std::pair<int,DecayTableEntry>(iPID,DecayTableEntry(iPID,im,ichainDecay)) );
    }
    channel randomDecay(int pID)
    {
        return table[pID].randomDecay();
    }     
    DecayTableEntry& operator[](int i){ return table[i];}  
};

// The main decay chain class. Each link (particle) in the decay chain is an instance of this class, with pointers to its mother and daughter links.
class chainParticle
{
    public:
        const double m; // Rest mass    
        // Constructor for the base node. Case: Base node is not a particle (but e.g. an annihilation).
        chainParticle(double E_COM, DecayTable *dc) : m(E_COM)
        {
            parent=NULL;
            child1=NULL;
            child2=NULL;
            p_parent=Ep4vec(vec3(0),E_COM); 
            pID=0;
            chainGeneration=0;
            decayTable = dc;
            isParent = false;
            boostToParentFrame = mat4::identity();
            boostToLabFrame    = mat4::identity();            
        }
        // Constructor for the base node. Case: Base node is a decaying particle.
        chainParticle(vec3 ipLab, DecayTable *dc, int ipID) : m((*dc)[ipID].m)
        {
            parent=NULL;
            child1=NULL;
            child2=NULL;
            p_parent=Ep4vec(ipLab,m);
            pID=ipID;
            chainGeneration=0;
            decayTable = dc;
            isParent = false;
            boostMatrixParentFrame(boostToParentFrame,p_parent);
            boostToLabFrame = boostToParentFrame;
        }   
        // Destructor
        ~chainParticle()
        {
            if(child1!=NULL)
                delete child1;
            if(child2!=NULL)
                delete child2;
        }
        // Iteratively add random links to the decay chain by Monte Carlo to a maximum length of maxSteps
        void generateDecayChainMC(int maxSteps)
        {
            if(!isParent && (*decayTable)[pID].chainDecay && chainGeneration < maxSteps)
            {
                channel chn = decayTable->randomDecay(pID);
                double m1 = (*decayTable)[chn.pID0].m;
                double m2 = (*decayTable)[chn.pID1].m;            
                double Etot = (*decayTable)[pID].m;
                double E1 = 0.5*(Etot*Etot+m1*m1-m2*m2)/Etot;
                double E2 = Etot-E1; 
                double abs_p = sqrt(E1*E1-m1*m1);                     
                vec3 dir = randOnSphere();
                vec4 p1(E1, abs_p*dir);
                vec4 p2(E2,-abs_p*dir);
                child1 = new chainParticle(p1, m1, decayTable, this, chainGeneration+1, chn.pID0);
                child2 = new chainParticle(p2, m2, decayTable, this, chainGeneration+1, chn.pID1);
                child1->generateDecayChainMC(maxSteps);
                child2->generateDecayChainMC(maxSteps);            
            }
            isParent = true;
        }
        // Add a new link to the decay chain by explicitly setting a decay for this particle
        void addChainLink(const channel &chn)
        {
            if(isParent)
            {
                cout << "Warning: Overwriting existing decay in decay chain." << endl;
                cutChain();
            }
            if(!(*decayTable)[pID].chainDecay || !(*decayTable)[pID].allowedDecay(chn))
            {
                cout << "Warning: Adding decay to decay chain that does not exist in decay table." << endl;
            }
            double m1 = (*decayTable)[chn.pID0].m;
            double m2 = (*decayTable)[chn.pID1].m;            
            double Etot = (*decayTable)[pID].m;
            double E1 = 0.5*(Etot*Etot+m1*m1-m2*m2)/Etot;
            double E2 = Etot-E1; 
            double abs_p = sqrt(E1*E1-m1*m1);                    
            vec3 dir = randOnSphere();
            const vec4 p1(E1, abs_p*dir);
            vec4 p2(E2,-abs_p*dir);
            child1 = new chainParticle(p1, m1, decayTable, this, chainGeneration+1, chn.pID0);
            child2 = new chainParticle(p2, m2, decayTable, this, chainGeneration+1, chn.pID1);           
            isParent = true;
        }    
        // Draw new angles for the decay products in this and all subsequent links of the decay chain
        void reDrawAngles()
        {
            if(isParent)
            {
                double m1 = child1->m;
                double m2 = child2->m;                
                vec3  dir = randOnSphere(); 
                double E1 = child1->p_parent[0];
                double E2 = child2->p_parent[0];    
                double abs_p = sqrt(E1*E1-m1*m1);                       
                vec4 p1(E1, abs_p*dir);
                vec4 p2(E2,-abs_p*dir);  
                child1->update(p1);
                child2->update(p2);                        
                child1->reDrawAngles();
                child2->reDrawAngles();  
            }
        }
        // Remove all subsequent links in the decay chain
        void cutChain()
        {
            if(isParent)
            {
                delete child1;
                delete child2;
                isParent = false;
            }
        }
        // Boost a given 4-momentum from this frame to the lab frame.
        vec4 p_to_Lab(const vec4 &p) const
        {
            return boostToLabFrame*p;
        }
        // Calculate 4-momentum of this particle in the lab frame.
        vec4 p_Lab() const
        {
            return parent->boostToLabFrame*p_parent;        
        }     
        // Iteratively collect all final state particles
        void collectFinalStates(vector<chainParticle*> &finalStates)
        {
            if(isParent)
            {
                child1->collectFinalStates(finalStates);
                child2->collectFinalStates(finalStates);                
            }       
            else 
            {
                finalStates.push_back(this);
            }
        }
        // Iteratively collect all final state particles of type ID
        void collectFinalStates(int ID, vector<chainParticle*> &finalStates)
        {
            if(pID == ID)
            {
                finalStates.push_back(this);
            }
            else
            {
                if(isParent)
                {
                    child1->collectFinalStates(ID,finalStates);
                    child2->collectFinalStates(ID,finalStates);                
                }
            }
        }    
        // Bracket operator. Elements 0 and 1 correspond to the two decay products.
        chainParticle& operator[](unsigned int i)
        {
            if(isParent && i<2)
                return (i==0) ? *child1 : *child2;
            else
            {
                cout << "Error: Trying to access non-existing decay chain entry." << endl;
                exit(1);
                return *this;
            }            
        }        
    private:  
        // Pointer to decay table
        DecayTable *decayTable;        
        // Particle properties
        mat4 boostToParentFrame;
        mat4 boostToLabFrame;
        vec4 p_parent;              // 4-momentum in parent's rest frame
        unsigned int pID;           // Particle identifier
        // Chain properties
        int chainGeneration;
        chainParticle *child1;
        chainParticle *child2;
        chainParticle *parent;    
        bool isParent;
        // Function for updating the Lorentz boost matrices according to a new 4-momentum.
        void update(vec4 &ip_parent)
        {
            p_parent = ip_parent;
            boostMatrixParentFrame(boostToParentFrame,p_parent);
            boostToLabFrame = parent->boostToLabFrame*boostToParentFrame;                 
        }
        // Constructor used by member functions during chain generation.
        chainParticle(const vec4 &pp, double im, DecayTable *dc, chainParticle *iparent, int ichainGeneration, int ipID) : m(im)
        {
            parent=iparent;
            child1=NULL;
            child2=NULL;
            p_parent=pp, 
            pID=ipID;
            chainGeneration=ichainGeneration;
            decayTable = dc;
            isParent = false;
            boostMatrixParentFrame(boostToParentFrame,p_parent);
            boostToLabFrame = parent->boostToLabFrame*boostToParentFrame;            
        }            
};

// Construct an energy-momentum 4-vector from the given 3-momentum and mass
vec4 Ep4vec(const vec3 p, double m)
{
    double E = sqrt(dot(p,p)+m*m);
    return vec4(E,p);
}

// Calculate the invariant mass of the given 4-vector pair
double invariantMass(const vec4 &a, const vec4 &b)
{
    vec4 tmp = a+b;
    return sqrt(dot(tmp,tmp));
}

typedef chainParticle decayChain;

// Test program
int main() 
{

    DecayTable dt;
    dt.addEntry(1, 1000, true);
    dt.addEntry(2, 500,  true);
    dt.addEntry(3, 0,   false);
    dt.addEntry(4, 100, false);   
    channel d1_23(2, 3, 0.5);
    channel d2_34(3, 4, 0.5); 
    dt[1].addChannel(d1_23);
    dt[2].addChannel(d2_34);   

    vec3 p1(5,5,5);
    decayChain chain(p1, &dt, 1);
    chain.addChainLink(d1_23);
    chain[0].addChainLink(d2_34);
    

    ofstream out("data", ios::out);
    if(!out.fail()) 
    {   
        for(int i=0;i<1000000;i++)
        {
            vec4 v1 = chain[1].p_Lab();
            vec4 v2 = chain[0][0].p_Lab();        
            chain.reDrawAngles();
            out << invariantMass(v1,v2) << endl;     
        }
    }
    out.close();       

    
    //double m1 = 1;
    //double m2 = 2;
    //double m3 = 41;
    //int pID = 2;
    //vec4 p1 = Ep4vec(vec3(5,5,5),m1);
    //vec4 p2 = Ep4vec(vec3(3,1,0),m2); 
    //vec4 p3 = Ep4vec(vec3(-5,71,-1),m3);
    
    
    //cout << "Invariant mass, frame 1: " << invariantMass(p1,p2) << endl;
    //cout << "Invariant mass, frame 2: " << invariantMass(p_parentFrame(p1, p3),p_parentFrame(p2, p3)) << endl;    
    //chainParticle x(p,m,pID);
    //vec4 pLab = x.p_Lab();
    //cout << p << endl;
    //cout << pLab << endl;
 
}
