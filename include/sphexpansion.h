/*
definitions for the SphExpansion class

MSP 24 Apr 2020 restructured
MSP  7 Oct 2020 fix density normalisation
MSP  4 Feb 2021 add dipole and quadrupole capability
MSP 25 May 2021 add ltrunc option; add self-gravitating coefficient calculation
MSP 28 Sep 2021 adjust radial order (nmax) truncation
MSP  9 Apr 2022 converted to Eigen
MSP 17 Jan 2023 expose coefficient tables for editing

*/

#ifndef SPHEXPANSION_H
#define SPHEXPANSION_H

#include <Eigen/Dense>
using Eigen::MatrixXd;

// MSP headers
// important preprocessor flags
#include "flags.h"

// converters from r to xi (for mapping tables)
#include "common/scaling.h"

// some basic basis elements
#include "common/basis.h"

// basic transformations from cartesian to spherical (and back)
#include "common/transform.h"

// basic translations from virial to physical units (and back)
#include "common/translate.h"

// the orientation stuff for centering the expansions
#include "common/orient.h"

// the cachefile stuff, also brings in the modelfile stuff
#include "sphere/sphcache.h"

// the coefficient stuff
#include "sphere/sphcoefs.h"

// set the default spherical harmonic evaluation flag
int SPHHARMONICDEFAULT = 2047;

#ifndef REPSVAL
double REPS=1.e-7;
#define REPSVAL
#endif


//using namespace std;
using std::cout, std::cerr, std::endl, std::setw, std::vector, std::ifstream, std::ios, std::string, std::ofstream, std::istringstream, std::make_tuple;

class SphExpansion
{
private:

  // doesn't need to be exposed to the outside world
  SphModel modeltable;
  SphCache cachetable; // does this actually have to be exposed?

  void initialise(string sph_cache_name,
      string model_file,
      string coef_file,
      string orient_file);

public:
  // the constructor
  SphExpansion(string sph_cache_name,
         string model_file,
         string coef_file,
         string orient_file="");

  // expose the coefficients for editing
  SphCoefs coeftable;

  // expose the important expansion data
  SphOrient orient;

  string sph_cache_name;
  string model_file;
  string coef_file;
  string orient_file;

  // expose the basic parameters of the expansion
  int LMAX; // maximum azimuthal order in expansion
  int NMAX; // maximum radial order in expansion
  int NUMT; // number of timesteps in the coefficients
  vector<double> T; // the array of timesteps from the coefficients
  // can use these for checking that different files match!

  // get potential function weights
  void get_pot_coefs(int l, int indx, int nmax, MatrixXd& coefs, MatrixXd& potd, MatrixXd& dpot, double *p, double *dp);

  // get density function weights
  void get_dens_coefs(int l, int indx, int nmax, MatrixXd& coefs, MatrixXd& dend, double *dd);


  // the base spherical class
  // harmonic flag is a bit flag for l orders (set as 2^l). see details in common/basis.h
  void determine_fields_at_point_sph(MatrixXd& coefs,
             double r, double theta, double phi,
             double& potl0, double& potl,
             double& potr, double& pott,
             double& potp,
             int harmonicflag=SPHHARMONICDEFAULT);

  // version with density return
  void determine_fields_at_point_sph(MatrixXd& coefs,
             double r, double theta, double phi,
             double& dens0, double& dens,
             double& potl0, double& potl,
             double& potr, double& pott,
             double& potp,
             int harmonicflag=SPHHARMONICDEFAULT);

  // version that is only density return
  void determine_fields_at_point_sph(MatrixXd& coefs,
             double r, double theta, double phi,
             double& dens0, double& dens,
             int harmonicflag=SPHHARMONICDEFAULT);

  // return the function weights at a given point
  std::tuple<MatrixXd,MatrixXd,MatrixXd,MatrixXd> determine_weights_at_point_sph(double r, double theta, double phi);

  // cartesian forces wrapper function
  void return_forces(MatrixXd& coefs,
         double x, double y, double z,
         double& fx, double& fy, double& fz,
         int harmonicflag=SPHHARMONICDEFAULT);

  // cartesian density wrapper function
  void return_density(MatrixXd& coefs,
          double x, double y, double z,
          double& d,
          int harmonicflag=SPHHARMONICDEFAULT);

  // coefficient interpolator
  //  this call may also be used as a coarse method to truncate coefficient series.
  //  however, it will not result in any speedup, as the coefficients are still evaluated (as zeros).
  void select_coefficient_time(double desired_time,
             MatrixXd& coefs_at_time, int ntrunc=-1, int ltrunc=0);


  // reset coefficients to the read-in values
  void reset_coefficients();

  std::tuple<std::vector<double>,std::vector<MatrixXd>> return_coefficients();

  void install_coefficients(std::vector<MatrixXd> tableau);

};

/*
std::tuple<vector<double>,std::vector<MatrixXd>> SphExpansion::output_coefficients()
{


  return make_tuple()
}
*/




SphExpansion::SphExpansion(string sph_cache_name,
         string model_file,
         string coef_file,
         string orient_file)
{
  SphExpansion::sph_cache_name = sph_cache_name;
  SphExpansion::model_file     = model_file;
  SphExpansion::coef_file      = coef_file;
  SphExpansion::orient_file    = orient_file;

  initialise(sph_cache_name, model_file, coef_file, orient_file);
}

void SphExpansion::initialise(string sph_cache_name,
                              string model_file,
                              string coef_file,
                              string orient_file="")
{
  // pull in the parts for the expansion
  try {
    read_sph_cache(sph_cache_name, SphExpansion::cachetable);
  } catch (const char* msg) {
    cerr << msg << endl;
    exit(1);
  }

  try {
    read_model(model_file, SphExpansion::modeltable);
  } catch (const char* msg) {
    cerr << msg << endl;
    exit(1);
  }

  try {
    read_coef_file (coef_file, SphExpansion::coeftable);
  } catch (const char* msg) {
    cerr << msg << endl;
    exit(1);
  }


  LMAX = coeftable.LMAX;
  NMAX = coeftable.NMAX;
  NUMT = coeftable.NUMT;

  // retrieve the internal time stamps (can this assignment fail?)
  T.resize(NUMT);
  T    = coeftable.t;

  // if no orient file, assume zeros?
  read_orient (orient_file, SphExpansion::orient);

  // finish setting up the model
  init_table(SphExpansion::modeltable, SphExpansion::cachetable);

}

void SphExpansion::reset_coefficients()
{
  /*
  reset coefficients back to those read in from the file
  */
  try {
    read_coef_file (SphExpansion::coef_file, SphExpansion::coeftable);
  } catch (const char* msg) {
    cerr << msg << endl;
    exit(1);
  }
}

std::tuple<std::vector<double>,std::vector<MatrixXd>> SphExpansion::return_coefficients()
{
  /*
  reset coefficients back to those read in from the file
  */
  return make_tuple(coeftable.t,coeftable.coefs);
}

void SphExpansion::install_coefficients(std::vector<MatrixXd> tableau)
{
  /*
  install external coefficients (presumably after some operations)
  */
  coeftable.coefs = tableau;
}


void SphExpansion::get_pot_coefs(int l, int indx, int nmax, MatrixXd& coefs, MatrixXd& potd, MatrixXd& dpot, double *p, double *dp)
{
  /*
    int l    : the harmonic order
    int indx : the indexed harmonic order (e.g. l and m values)

   */
  double pp, dpp;
  //int i;

  pp = dpp = 0.0;

  for (int i=0; i<nmax; i++) {
    pp  += potd(l,i) * coefs(indx,i);
    dpp += dpot(l,i) * coefs(indx,i);
  }

  *p = -pp;
  *dp = -dpp;
}

void SphExpansion::get_dens_coefs(int l, int indx, int nmax, MatrixXd& coefs, MatrixXd& dend, double *dd)
{
  /*
    int l    : the harmonic order
    int indx : the indexed harmonic order (e.g. l and m values)

   */
  double daccum;
  //int i;

  daccum = 0.0;

  for (int i=0; i<nmax; i++)
    daccum  += dend(l,i) * coefs(indx,i);

  *dd = daccum;
}


//void SphExpansion::determine_weights_at_point_sph(MatrixXd& coefs,
//                                                  double r, double theta, double phi,
//                                                  MatrixXd& potl,
//                                                  MatrixXd& potr, MatrixXd& pott, MatrixXd& potp)
std::tuple<MatrixXd,MatrixXd,MatrixXd,MatrixXd> SphExpansion::determine_weights_at_point_sph(double r, double theta, double phi)
{
   /*
   return just the function weights at a position

   */

   MatrixXd potl,potr,pott,potp;
   potl.resize((LMAX+1)*(LMAX+1),NMAX);
   potr.resize((LMAX+1)*(LMAX+1),NMAX);
   pott.resize((LMAX+1)*(LMAX+1),NMAX);
   potp.resize((LMAX+1)*(LMAX+1),NMAX);

   int numl = cachetable.LMAX;

   int l,loffset,moffset,m,n;
   double rs,fac1,fac2,fac3,fac4,costh,dp;
   double p,pc,dpc,ps,dps;

   // block here, some problem with a zero in theta here. TBD.
   if (theta<1.e-6) theta = 1.e-6;
   costh = cos(theta);

   // block zero evaluation
   if (r<REPS) r = REPS;

   MatrixXd potd,dpot;
   get_dpotl(r, cachetable, potd, dpot);

   MatrixXd factrl;
   factorial(numl, factrl);

   MatrixXd legs, dlegs;
   dlegendre_R(numl, costh, legs, dlegs);

   vector<double> cosm(cachetable.NMAX),sinm(cachetable.NMAX);
   sinecosine_R(numl, phi, cosm, sinm);

   for (l=0, loffset=0; l<=numl; loffset+=(2*l+1), l++) {

     // m loop
     for (m=0, moffset=0; m<=l; m++) {

       if (l==0) { fac1 = 0.25/M_PI; } else { fac1 = (2.0*l+1.0)/(4.0*M_PI);};

       if (m==0) {
         fac2 = fac1*legs(l,m);

         for (n=0; n<cachetable.NMAX; n++) {
           potl(loffset+moffset,n) = fac2*potd(l,n);
           potr(loffset+moffset,n) = fac2*dpot(l,n);
           pott(loffset+moffset,n) = fac1*dlegs(l,m)*potd(l,n);
           potp(loffset+moffset,n) = 0.0;
         }

         moffset++;

       }
       else {
         fac2 = 2.0 * fac1 * factrl(l,m);
         fac3 = fac2 *  legs(l,m);
         fac4 = fac2 * dlegs(l,m);

         for (n=0; n<cachetable.NMAX; n++) {
           potl(loffset+moffset,n) = fac3*( potd(l,n)*cosm[m]);
           potr(loffset+moffset,n) = fac3*( dpot(l,n)*cosm[m]);
           pott(loffset+moffset,n) = fac4*( potd(l,n)*cosm[m]);
           potp(loffset+moffset,n) = fac3*(-potd(l,n)*sinm[m])*m;
           potl(loffset+moffset+1,n) = fac3*( potd(l,n)*sinm[m]);
           potr(loffset+moffset+1,n) = fac3*( dpot(l,n)*sinm[m]);
           pott(loffset+moffset+1,n) = fac4*( potd(l,n)*sinm[m]);
           potp(loffset+moffset+1,n) = fac3*( potd(l,n)*cosm[m])*m;
         }

         moffset +=2;
       }
     }
   }
   potl *= -1;
   potr *= -1;
   pott *= -1;
   potp *= -1;

   return make_tuple(potl,potr,pott,potp);
}


void SphExpansion::determine_fields_at_point_sph(MatrixXd& coefs,
                                                 double r, double theta, double phi,
                                                 double& potl0, double& potl,
                                                 double& potr, double& pott, double& potp,
                                                 int harmonicflag)
{
  /*
  version without density

  @IMPROVE: remove potl0

  see the equivalent exp call in SphericalBasis.cc

  */

  int numl = cachetable.LMAX;

  int l,loffset,moffset,m;
  double rs,fac1,fac2,fac3,fac4,costh,dp;
  double p,pc,dpc,ps,dps;

  // block here, some problem with a zero in theta here. TBD.
  if (theta<1.e-6) theta = 1.e-6;
  costh = cos(theta);

  // block zero evaluation
  if (r<REPS) r = REPS;

  fac1 = 0.25/M_PI;

  MatrixXd potd,dpot;
  get_dpotl(r, cachetable, potd, dpot);

  // is this ever evaluating the l=0,n>0 terms??

  // compute the monopole values
  get_pot_coefs(0, 0, cachetable.NMAX, coefs, potd, dpot, &p, &dp);
  potl = fac1*p;
  potr = fac1*dp;
  pott = potp = 0.0;

  // l loop
  if (harmonicflag==0) return;

  MatrixXd factrl;
  factorial(numl, factrl);

  MatrixXd legs, dlegs;
  dlegendre_R(numl, costh, legs, dlegs);

  vector<double> cosm(cachetable.NMAX),sinm(cachetable.NMAX);
  sinecosine_R(numl, phi, cosm, sinm);

  for (l=1, loffset=1; l<=numl; loffset+=(2*l+1), l++) {

    if (check_flags(harmonicflag,l)==0) continue;


    // m loop
    for (m=0, moffset=0; m<=l; m++) {
      fac1 = (2.0*l+1.0)/(4.0*M_PI);
      if (m==0) {
        fac2 = fac1*legs(l,m);

        get_pot_coefs(l, loffset+moffset, cachetable.NMAX, coefs, potd, dpot, &p, &dp);
        potl += fac2*p;
        potr += fac2*dp;
        pott += fac1*dlegs(l,m)*p;
        moffset++;
      }
      else {
        fac2 = 2.0 * fac1 * factrl(l,m);
        fac3 = fac2 *  legs(l,m);
        fac4 = fac2 * dlegs(l,m);

        get_pot_coefs(l, loffset+moffset,   cachetable.NMAX, coefs, potd, dpot, &pc, &dpc);
        get_pot_coefs(l, loffset+moffset+1, cachetable.NMAX, coefs, potd, dpot, &ps, &dps);

        potl += fac3*( pc*cosm[m] + ps*sinm[m]);
        potr += fac3*(dpc*cosm[m] + dps*sinm[m]);
        pott += fac4*( pc*cosm[m] +  ps*sinm[m]);
        potp += fac3*(-pc*sinm[m] +  ps*cosm[m])*m;

        moffset +=2;
      }
    }
  }

  // zero out pott/potp if below resolution limit for forces
  // (e.g. centre crossing problems)
  if (r<cachetable.RMIN) {
    pott = 0.;
    potp = 0.;
  }


}


void SphExpansion::determine_fields_at_point_sph(MatrixXd& coefs,
 double r, double theta, double phi,
 double& dens0, double& dens,
 double& potl0, double& potl,
 double& potr, double& pott, double& potp,
 int harmonicflag)
{
  /*
  // version WITH density

  // no potl0/dens0 definition?

  see the equivalent exp call in SphericalBasis.cc

  */


  int numl = cachetable.LMAX;

  int l,loffset,moffset,m;
  double rs,fac1,fac2,fac3,fac4,costh,dp;
  double p,pc,dpc,ps,dps,d,dc,ds;

  // block here, some problem with a zero in theta here. TBD.
  if (theta<1.e-6) theta = 1.e-6;
  costh = cos(theta);

  // block zero evaluation
  if (r<REPS) r = REPS;

  //fac0 = 4.*M_PI;
  fac1 = 0.25/M_PI;

  MatrixXd potd,dpot,dend;
  get_dpotl_density(r, cachetable, potd, dpot, dend);

  // is this ever evaluating the l=0,n>0 terms??

  // compute the monopole values
  get_pot_coefs(0, 0, cachetable.NMAX, coefs, potd, dpot, &p, &dp);
  potl = fac1*p;
  potr = fac1*dp;
  pott = potp = 0.0;

  get_dens_coefs(0, 0, cachetable.NMAX, coefs, dend, &d);
  dens  = fac1*d;
  dens0 = fac1*d;


  // l loop
  if (harmonicflag==0) return;

   MatrixXd factrl;
  factorial(numl, factrl);

  MatrixXd legs, dlegs;
  dlegendre_R(numl, costh, legs, dlegs);

  vector<double> cosm(cachetable.NMAX),sinm(cachetable.NMAX);
  sinecosine_R(numl, phi, cosm, sinm);

  for (l=1, loffset=1; l<=numl; loffset+=(2*l+1), l++) {

    if (check_flags(harmonicflag,l)==0) continue;

    // m loop
    for (m=0, moffset=0; m<=l; m++) {
      fac1 = (2.0*l+1.0)/(4.0*M_PI);
      if (m==0) {
  fac2 = fac1*legs(l,m);

  get_dens_coefs(l,loffset+moffset,cachetable.NMAX, coefs, dend, &d);
  dens += fac2*d;

  get_pot_coefs(l, loffset+moffset, cachetable.NMAX, coefs, potd, dpot, &p, &dp);
  potl += fac2*p;
  potr += fac2*dp;
  pott += fac1*dlegs(l,m)*p;
  moffset++;
      }
      else {
  fac2 = 2.0 * fac1 * factrl(l,m);
  fac3 = fac2 *  legs(l,m);
  fac4 = fac2 * dlegs(l,m);

  get_dens_coefs(l,loffset+moffset,  cachetable.NMAX, coefs, dend, &dc);
  get_dens_coefs(l,loffset+moffset+1,cachetable.NMAX, coefs, dend, &ds);
  dens += fac3*(dc*cosm[m] + ds*sinm[m]);

  get_pot_coefs(l, loffset+moffset,   cachetable.NMAX, coefs, potd, dpot, &pc, &dpc);
  get_pot_coefs(l, loffset+moffset+1, cachetable.NMAX, coefs, potd, dpot, &ps, &dps);

  potl += fac3*( pc*cosm[m] + ps*sinm[m]);
  potr += fac3*(dpc*cosm[m] + dps*sinm[m]);
  pott += fac4*( pc*cosm[m] +  ps*sinm[m]);
  potp += fac3*(-pc*sinm[m] +  ps*cosm[m])*m;

  moffset +=2;
      }
    }
  }

  // zero out pott/potp if below resolution limit for forces
  // (e.g. centre crossing problems)
  if (r<cachetable.RMIN) {
    pott = 0.;
    potp = 0.;
    dens = cachetable.d0[0]; // set to the smallest value of density. check if this is a good idea
  }


}



void SphExpansion::determine_fields_at_point_sph(MatrixXd& coefs,
 double r, double theta, double phi,
 double& dens0, double& dens,
 int harmonicflag)
{
  /*
  // version that is ONLY density

  see the equivalent exp call in SphericalBasis.cc

  */

  int numl = cachetable.LMAX;

  int l,loffset,moffset,m;
  double rs,fac1,fac2,fac3,fac4,costh,dp;
  double p,pc,dpc,ps,dps,d,dc,ds;

  // block here, some problem with a zero in theta here. TBD.
  if (theta<1.e-6) theta = 1.e-6;
  costh = cos(theta);

  // block zero evaluation
  if (r<REPS) r = REPS;

  fac1 = 0.25/M_PI;

  // is this ever evaluating the l=0,n>0 terms??
  MatrixXd potd,dpot,dend;
  get_dpotl_density(r, cachetable, potd, dpot, dend);

  // compute the monopole values
  get_dens_coefs(0, 0, cachetable.NMAX, coefs, dend, &d);
  dens  = fac1*d;
  dens0 = fac1*d;

  // l loop
  if (harmonicflag==0) return;

   MatrixXd factrl;
  factorial(numl, factrl);

  MatrixXd legs, dlegs;
  dlegendre_R(numl, costh, legs, dlegs);

  vector<double> cosm(cachetable.NMAX),sinm(cachetable.NMAX);
  sinecosine_R(numl, phi, cosm, sinm);

  for (l=1, loffset=1; l<=numl; loffset+=(2*l+1), l++) {

    if (check_flags(harmonicflag,l)==0) continue;

    // m loop
    for (m=0, moffset=0; m<=l; m++) {
      fac1 = (2.0*l+1.0)/(4.0*M_PI);
      if (m==0) {
        fac2 = fac1*legs(l,m);

        get_dens_coefs(l,loffset+moffset,cachetable.NMAX, coefs, dend, &d);
        dens += fac2*d;

        moffset++;
      }
      else
      {
        fac2 = 2.0 * fac1 * factrl(l,m);
        fac3 = fac2 *  legs(l,m);

        get_dens_coefs(l,loffset+moffset,  cachetable.NMAX, coefs, dend, &dc);
        get_dens_coefs(l,loffset+moffset+1,cachetable.NMAX, coefs, dend, &ds);
        dens += fac3*(dc*cosm[m] + ds*sinm[m]);

        moffset +=2;
      }
    }
  }

  // zero out pott/potp if below resolution limit for forces
  // (e.g. centre crossing problems)
  if (r<cachetable.RMIN) {
    dens = cachetable.d0[0]; // set to the smallest value of density. check if this is a good idea
  }


}




void SphExpansion::return_forces(MatrixXd& coefs,
                                 double x, double y, double z,
                                 double& fx, double& fy, double& fz,
                                 int harmonicflag)
{
  /*
    test force return from just one component, from the centre of the expansion
   */

  // translate all times and positions into exp virial units
  double xvir,yvir,zvir;
  physical_to_virial_length(x,y,z,xvir,yvir,zvir);

  double rtmp,phitmp,thetatmp;
  double tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  cartesian_to_spherical(xvir, yvir, zvir, rtmp, phitmp, thetatmp);

  determine_fields_at_point_sph(coefs, rtmp,thetatmp,phitmp,
                                tpotl0,tpotl,
                                fr,ft,fp,harmonicflag);

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fx,fy,fz);

}


void SphExpansion::return_density(MatrixXd& coefs,
                                  double x, double y, double z,
                                  double& d,
                                  int harmonicflag)
{
  /*
    return density
   */

  // translate all times and positions into exp virial units
  double xvir,yvir,zvir;
  physical_to_virial_length(x,y,z,xvir,yvir,zvir);

  double rtmp,phitmp,thetatmp;
  double tpotl0,tpotl,fr,ft,fp,tdens0;
  d = 0;

  cartesian_to_spherical(xvir, yvir, zvir, rtmp, phitmp, thetatmp);

  determine_fields_at_point_sph(coefs,
        rtmp,thetatmp,phitmp,
        tdens0,d,harmonicflag);

#if DEEPDEBUGCOEFS
  std::cout << "sphexpansion.h::return_density:"  << setw(14) << rtmp << setw(14) << thetatmp << setw(14) << phitmp << setw(14) << fr << setw(14) << ft << setw(14) << fp << '\n';
#endif

}

void SphExpansion::select_coefficient_time(double desired_time,
           MatrixXd& coefs_at_time, int ntrunc, int ltrunc) {
  /*
    linear interpolation to get the coefficient matrix at a specific time

   time units must be virial time units to match the input coefficient
   table

   if ntrunc<0, will select all orders. if ntrunc>=0, will truncate at
   the specified n_max by setting all coefficients above order to 0.0.

   the ntrunc limit will be applied for all orders at or above ltrunc.
   by default, the ntruncation will be applied for all l orders.
   */

  int numl;

  numl = (coeftable.LMAX+1)*(coeftable.LMAX+1);

  coefs_at_time.resize(numl,coeftable.NMAX);



  // starting at the first indx, stop when we get to the matching time
  int indx = 0;
  while (coeftable.t[indx]<=desired_time) {
    indx ++;
  }

  // check the spacing on coeftable.t (can be nonuniform)
  double dt = coeftable.t[indx+1] - coeftable.t[indx];

  // guard against wanton extrapolation: should this stop the model?
  //if (indx<0) cerr << "select_coefficient_time: time prior to simulation start selected. setting to earliest step." << endl;

  // guard against going past the end of the simulation
  if (indx>coeftable.NUMT-2) {
    cerr << "select_coefficient_time: time after to simulation end selected. setting to latest step, but this behaviour is not supported." << "\n";
    indx = coeftable.NUMT-2;
  }

#if DEEPDEBUGTIME
  std::cout << "indx=" << indx
       << " MinT=" << coeftable.t[0]
       << " MaxT=" << coeftable.t[coeftable.NUMT-1]
       << " desired_time=" << desired_time
       << std::endl;
#endif

  if (indx<=0) {

#if DEEPDEBUGTIME
  std::cout << "sphexpansion.h: starting before the simulation starts..."
       << std::endl;
#endif
    // set the coefficients to be the first coefficients from the simulation

    // do these need to be loops now that we have eigen?
    coefs_at_time = coeftable.coefs[0];

    /*
    for (int l=0; l<numl; l++){
      for (int n=0; n<coeftable.NMAX; n++) {
        coefs_at_time(l,n) = coeftable.coefs[0](l,n);
      }
    }
    */

  } else {

    // case where the simulations are not before the beginning of the simulation
    // interpolate from the two closest times

    double x1 = (coeftable.t[indx+1] - desired_time)/dt;
    double x2 = (desired_time - coeftable.t[indx])/dt;

#if DEEPDEBUGTIME
    cout << "dt=" << setw(16) << dt << "  t[indx+1]="  << setw(16) << coeftable.t[indx+1]
                              << "  t[indx]="  << setw(16) << coeftable.t[indx]
                                    << "  desiredT="  << setw(16) << desired_time << "\n";
    cout << "x1=" << setw(16) << x1 << " x2=" << setw(14) << x2 << "\n";
#endif

    // do some eigen multiplication to get the coefficient table. this is faster, right?
    coefs_at_time = (x1 * coeftable.coefs[indx] + x2 * coeftable.coefs[indx+1]);
    /*
    for (int l=0; l<numl; l++){
      for (int n=0; n<coeftable.NMAX; n++) {
        coefs_at_time(l,n) = (x1 * coeftable.coefs[indx](l,n) + x2 * coeftable.coefs[indx+1](l,n));
      }
    }
    */
  }

  // go through and zero out non-selected orders
  //   note: will be applied for all l orders >= ltrunc

  if (ntrunc>=0) {
    for (int l=ltrunc*ltrunc; l<numl; l++){

      // start at order+1 so this call will work if n=0.
      for (int n=ntrunc+1; n<coeftable.NMAX; n++) {
        coefs_at_time(l,n) = 0.0;
      }

    }
  }

}

#endif
