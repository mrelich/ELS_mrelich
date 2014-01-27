//===================================================================
//    Electron Beam Initial Information :: PrimaryEbeam.cc
//      Author    : T.Shibata
//      Creation  : 2006.02.13
//  Last Update   : 2010.02.11
//===================================================================
#include "PrimaryEbeam.hh"
//-----------------------------------------------
using namespace std;
//------------------------------------------------
Ebeam::Ebeam(){

  time_t  t;
  //srand(time(&t) % RAND_MAX);
  srand( (unsigned) time(NULL) );

  pi = 3.1415926535;

  Nppb=10;
  bunchtime=0.35;   // unit = nsec 0.35nsec = 350ps = S-band wavelenght

}
//------------------------------------------------
Ebeam::~Ebeam(){}
//------------------------------------------------
double Ebeam::rn(void)
{
  return double(rand())/double(RAND_MAX);
}
//------------------------------------------------
double Ebeam::gaussian( double mean, double sigma )
{
  double rn1 = rn();
  double rn2 = rn();
  double rrn = rn();

  double eta1 = sigma*sqrt(-2*log(rn1))*cos(2*pi*rn2);
  double eta2 = sigma*sqrt(-2*log(rn1))*sin(2*pi*rn2);

  double newx(0.0);
  if(rrn >= 0.5){
    newx=eta1;
  }else{
    newx=eta2;
  }
  return mean+newx;
}
//------------------------------------------------
double Ebeam::double_gaussian( double fluction,
		      	       double mean1, double sigma1,
			       double mean2, double sigma2 )
{
  return fluction*gaussian(mean1,sigma1)+(1.0-fluction)*gaussian(mean2,sigma2);
}
//------------------------------------------------
double Ebeam::rn_radian( void )
{  
  return rn()*2.0*pi;
}                                                      
//------------------------------------------------
double Ebeam::pulse( double xs, double xe )
{
  return xs + rn()*( xe -xs );
}
//------------------------------------------------
int Ebeam::poisson( double mean, double xmax, double xmin )
{
  double rnx, rny, p ;
  int irnx;

  while(1)
    {
      rnx = ( xmax - xmin )*rn();
      rny = rn();
      irnx = int( rnx );
      p = probability_poisson( irnx, mean );
      if( rny < p ) break;
    }
  return irnx;
}
//------------------------------------------------
double Ebeam::probability_poisson( int n, double m )
{
  double tmp1=1;
  double tmp2=1;
  for( int i=0 ; i<n; i++ ){
    tmp1*=m;
    tmp2*=(i+1);
  }
  return exp(-m)*tmp1/tmp2;
}
//------------------------------------------------
double Ebeam::bunch_time( int Nevent )
{
  int i = Nevent/Nppb;
  return double(i)*bunchtime;
}
//------------------------------------------------
double Ebeam::crystal_ball( double mean,
			    double sigma,
			    double alpha,
                            double n,
                            double xmin, double xmax )
{ 

  double x(0.);
  double P1(0.);
  double P2(0.);

  while(1){
        x  = rn()*( xmax - xmin ) + xmin;
        P1 = probability_crystal_ball( x, mean, sigma, alpha, n );
        P2 = rn();
	if( P2 < P1 ) break;
  }

  return x;
}
//------------------------------------------------
double Ebeam::probability_crystal_ball( double x,
                                        double mean,
                                        double sigma,
                                        double alpha,
                                        double n )
{
  double rt;

  double A = pow( n/fabs(alpha), n )*exp( (-1)*alpha*alpha/2.0 );
  double B = n/fabs(alpha) - fabs(alpha);

  if ( (x-mean)/sigma > (-1)*alpha ){
    rt = exp( (-1)*(x-mean)*(x-mean)/(2.0*sigma*sigma) );
  }else{
    rt = A*pow( B - (x-mean)/sigma , (-1)*n );
  }

  return rt; // return probability
}
//------------------------------------------------
double Ebeam::emittance_position( double emmitance[], double mean )
{
  double gamma = emmitance[0];
  double beta  = emmitance[1];
  double alpa  = emmitance[2];
  double em    = emmitance[3];
  gamma = (1.+alpa*alpa)/beta;

  double sigma = sqrt(beta*em);   // unit = m
  return pulse( mean-sigma, mean+sigma );
}
//------------------------------------------------
double Ebeam::emittance_direction( double emmitance[], double posi )
{
  double gamma = emmitance[0];
  double beta  = emmitance[1];
  double alpa  = emmitance[2];
  double em    = emmitance[3];
  gamma = (1.+alpa*alpa)/beta;

  double mean  = (-1)*alpa/beta*posi; // unit = m
  double sigma2 = mean*mean - (gamma*posi*posi-em)/beta;
  double sigma = sqrt( sigma2 );

  return pulse( mean-sigma, mean+sigma );
}
//------------------------------------------------

//------------------------------------------------
// Special Beam Direction : 2010.11.02
// For Study the Beam Energy Spectrum 
//------------------------------------------------
void Ebeam::special_beam_direction_20101102( double beam_position[],
	                   		     double &beam_direction_x,
					     double &beam_direction_y,
                                             double &beam_direction_z )
{
  double target_position[3];
 	   double shift_r = 5.0*rn();
           double shift_t = rn_radian();
           double shift_y = shift_r*cos(shift_t);   
	   double shift_z = shift_r*sin(shift_t);
           target_position[0]=11131.3;
  	   target_position[1]=-1718.2 + shift_y;
	   target_position[2]=1300.3  + shift_z;
     
  double beam_direction[3];
  for( int i=0; i<3; i++ ) beam_direction[i]=target_position[i]-beam_position[i];
  double mag=sqrt( beam_direction[0]*beam_direction[0]
                 + beam_direction[1]*beam_direction[1]
          	 + beam_direction[2]*beam_direction[2] );

  for( int i=0; i<3; i++ ) beam_direction[i]/=mag;
  beam_direction_x=beam_direction[0];
    beam_direction_y=beam_direction[1];
      beam_direction_z=beam_direction[2];

  return;
}

