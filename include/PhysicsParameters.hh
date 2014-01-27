//===================================================================
//    Physics Parameters
//      Author    : T.Shibata
//
//      Creation  : 2006.03.08
//    Last Update : 2012.09.24
//
//===================================================================

// Math  Value ///////////////////

static const double PI     ( 3.14159265358979 );
static const double PI2    ( 2.0*PI );

static const double RADI   ( 2.0*PI/360.0 );

// Physics Value ////////////////

//static const double  c_light (2.99792458E8 );           // m/s
static const double  Kelvin(-273.15);                  // K

static const double  hc    (197.327053);                // MeV * fm

static const double  qe    (1.602176462E-19 );          // C
static const double  k_boltzman( 1.3806503E-23   );     // m^2 * kg * s^-2 * K^-1 ( = J/K )
static const double  N_avogadro( 6.02214199E+23  );     // mol^-1

// 
static const double ft2mm ( 304.8 ); // ft to mm 
static const double inch2mm ( 25.4 ); // inch to mm 

// added in 2012.09.24
static double PhotonLambdatoEnergy( double lambda ){ return 2.0*PI*hc/lambda; }
static double PhotonEnergytoLambda( double energy ){ return 2.0*PI*hc/energy; }

