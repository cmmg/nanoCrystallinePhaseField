//mesh parameters
#define DIMS 2
#define refinementFactor 8
#define problemWidth 1.0

// mechanics and solte drag control
#define isMechanics true
#define isSoluteDrag true
#define isFiniteStrain true

//mechanics and solute drag increment
#define mechanicsStartIncrement 20
#define mechanicsEndIncrement 30
#define dragStartIncrement 40
#define dragEndIncrement 500

// grain growth parameters
#define N_seed_points 60 //(no of grain seed points)
#define n_diff_grains 6  //(no of order parameters)
#define InterfaceEnergyParameter 1.0e-3 //(interface energy  parameter for grain boundary)

// Solute parameters

#define WA 1.0
#define WB -0.1  
#define kappa 5.0e-3 
#define n_solute 1 //(degrees of freedom for solute particles per node)
#define n_chemical_potential 1 //(degrees of freedom for chemical potential per node)

//kinetic parameters                                                             
#define TimeStep 5.0e-7
#define TotalTime 11000*TimeStep
#define Mobility 100.0e2 //(order parameter mobility)
#define M_alpha 1.0 //(solute mobility)

//other parameters and variables
#define PI 3.1415
#define outputFileName "solution"
#define Vm 1.0
#define wellHeight 1.0

// elastic modulii
#define alpha1 2000 //(E1 along e1-direction)
#define beta1 1000  //(E2 along e2-direction)

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
