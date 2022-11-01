//mesh parameters
#define DIMS 3
#define refinementFactor 7
#define maxRefinementLevel 7
#define minRefinementLevel 3
#define problemWidth 1.0

// mechanics and solte drag control
#define isMechanics false
#define isSoluteDrag false
#define isFiniteStrain false

//mechanics and solute drag increment                                                                                
#define mechanicsStartIncrement 20
#define mechanicsEndIncrement 30
#define dragStartIncrement 40
#define dragEndIncrement 1000

// grain growth parameters
#define N_seed_points 60
#define n_diff_grains 6 //4
#define InterfaceEnergyParameter 1.0e-3

// Solute parameters
#define WA 1.0
#define WB 0.1
#define kappa 5.0e-3
#define n_solute 1
#define n_chemical_potential 1

//kinetic parameters
#define TimeStep 5.0e-7
#define TotalTime 11000*TimeStep
#define Mobility 100.0e2
#define M_alpha 0.5

//other parameters and variables
#define PI 3.1415
#define outputFileName "solution"
#define Vm 1.0
#define wellHeight 1.0

// elastic modulii
#define alpha1 2000
#define beta1 1000

// degree of freedom per node
#define TotalDOF n_diff_grains

#if isMechanics
#undef TotalDOF
#define TotalDOF n_diff_grains + DIMS
#endif

#if isSoluteDrag
#undef TotalDOF
#define TotalDOF n_diff_grains + n_solute + n_chemical_potential
#endif

#if (isMechanics && isSoluteDrag)
#undef TotalDOF
#define TotalDOF DIMS + n_diff_grains + n_solute + n_chemical_potential
#endif
