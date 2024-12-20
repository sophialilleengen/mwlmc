/*
header file for fullintegrate.cc capabilities

MSP 21 Apr 2022 cleaned version v0.2.3

*/
#ifndef FULLINTEGRATE_H
#define FULLINTEGRATE_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdio.h>

// eigen includes
#include <Eigen/StdVector>
#include <Eigen/Dense>
using Eigen::MatrixXd;

// set model parameters
//#include "modelfilesLilleengen2023.h"     // this file reproduces the Lilleengen+ 2023 results
//#include "modelfilesLilleengen2023m0.h"   // this file is a test version of the Lilleengen+ 2023 model, with m=0 restriction in the disc
//#include "modelfilesErkal2019m0.h"            // this file is a model re-scaled to try and match the Erkal+ 2019 model more closely

#include "modelfilesErkal2019.h"            // this file is a model re-scaled to try and match the Erkal+ 2019 model more closely


// the spherical expansion headers
#include "sphexpansion.h"

// the cylindrical expansion headers
#include "cylexpansion.h"

// small value
#ifndef EEPSVAL
double EEPS=1.e-4;
#define EEPSVAL
#endif


class MWLMC
{
private:

  void initialise();

  // version for integration, to avoid extra allocations
  void all_forces_coefs(MatrixXd mwcoefs, MatrixXd lmccoefs, MatrixXd mwdcoscoefs, MatrixXd mwdsincoefs,
                        double t, double x, double y, double z,
                        double& fx, double& fy, double& fz,
                        int mwhharmonicflag=127, int mwdharmonicflag=127, int lmcharmonicflag=127,
                        bool verbose=false);

  // if integrating just the MW
  void mw_forces_coefs(MatrixXd mwcoefs,
                       MatrixXd mwdcoscoefs, MatrixXd mwdsincoefs,
                       double t, double x, double y, double z,
                       double& fx, double& fy, double& fz,
                       int mwhharmonicflag, int mwdharmonicflag,
                       bool verbose=false);

  // get centres of the expansions with VIRIAL time input
  std::vector<double> get_expansion_centres_virial(double tvir, bool verbose=false);
  std::vector<double> get_expansion_centre_velocities_virial(double tvir, bool verbose=false);

  std::vector<double> get_lmc_centre_virial(double tvir, bool verbose);

  MatrixXd get_trajectories(double dt=native_timestep, bool virial=false);

public:

  // the expansion classes: these might be better private
  SphExpansion* MW;
  SphExpansion* LMC;
  CylExpansion* MWD;

  // the constructor (no type, no arguments, all defined in modelfiles.h)
  MWLMC();

  // return all fields for the MW halo (in the frame of the MW halo)
  std::vector<double> mwhalo_fields(double t, double x, double y, double z,
                                    //bool globalframe=false,
                                    int mwhharmonicflag=127, bool verbose=false);

  MatrixXd mwhalo_fields(double t, std::vector<double> x, std::vector<double> y, std::vector<double> z,
                         //bool globalframe=false,
                         int mwhharmonicflag=127, bool verbose=false);

  // return all fields for the LMC halo
  std::vector<double> lmc_fields(double t, double x, double y, double z,
                                 //bool globalframe=false,
                                 int lmcarmonicflag=127, bool verbose=false);

  MatrixXd lmc_fields(double t, std::vector<double> x, std::vector<double> y, std::vector<double> z,
                      //bool globalframe=false,
                      int lmcarmonicflag=127, bool verbose=false);

  // return all fields for the MW disc
  // NOTE: density does not work here. not enabled yet for cylindrical expansions.
  //       leave as an inspirational placeholder
  std::vector<double> mwd_fields(double t, double x, double y, double z,
                                 //bool globalframe=false,
                                 int mwdharmonicflag=127, bool verbose=false);

  MatrixXd mwd_fields(double t, std::vector<double> x, std::vector<double> y, std::vector<double> z,
                      //bool globalframe=false,
                      int mwdharmonicflag=127, bool verbose=false);

  // return total forces
  // @IMPROVE: write all potential as well
  std::vector<double> all_forces(double t, double x, double y, double z,
                                 //bool globalframe=true,
                                 int mwhharmonicflag=127, int mwdharmonicflag=127, int lmcharmonicflag=127,
                                 bool verbose=false);

  MatrixXd all_forces(double t, std::vector<double> x, std::vector<double> y, std::vector<double> z,
                      //bool globalframe=true,
                      int mwhharmonicflag=127, int mwdharmonicflag=127, int lmcharmonicflag=127,
                      bool verbose=false);

  // get centres of the expansions with PHYSICAL time input
  std::vector<double> get_expansion_centres_physical(double t, bool verbose=false);
  std::vector<double> get_expansion_centre_velocities_physical(double t, bool verbose=false);

  // compute an orbit integration in all three components
  // this uses PHYSICAL units by nature
  // lengths in kpc
  // velocities in km/s
  // time in Gyr

  // compute an orbit integration using only the initial Milky Way potential
  MatrixXd mworbit(vector<double> xinit,
                   vector<double> vinit,
                   double tbegin=-2.5,
                   double tend=0.0,
                   double dt=0.002);

   // compute an orbit integration using only the initial Milky Way potential
   std::vector< MatrixXd > mworbit(MatrixXd xinit,
                                   MatrixXd vinit,
                                   double tbegin=-2.5,
                                   double tend=0.0,
                                   double dt=0.002);

   // compute an orbit rewind in all three components
   MatrixXd rewind(vector<double> xinit,
                   vector<double> vinit,
                   double dt=0.002,
                   int mwhharmonicflag=127, int mwdharmonicflag=127, int lmcharmonicflag=127,
                   double rewind_time=2.5,
                   bool discframe = true,
                   bool verbose = false);

   std::vector< MatrixXd > rewind(MatrixXd xinit,
                                  MatrixXd vinit,
                                  double dt=0.002,
                                  int mwhharmonicflag=127, int mwdharmonicflag=127, int lmcharmonicflag=127,
                                  double rewind_time=2.5,
                                  bool discframe = true,
                                  bool verbose = false);

  // get the LMC trajectory (relative to the MW disc centre)
  MatrixXd get_lmc_trajectory(double rewindtime=2.5,double dt=native_timestep);

  // print an orbit array
  void print_orbit(MatrixXd orbit, string orbitfile);

  // reset coefficients to cache values
  void reset_mw_coefficients();
  void reset_all_coefficients();

  std::tuple<MatrixXd,MatrixXd,MatrixXd,MatrixXd> get_mw_function_weights(double x, double y, double z);
  std::tuple<std::vector<double>,std::vector<MatrixXd>> return_mw_coefficients();
  void install_mw_coefficients(std::vector<MatrixXd> tableau);

  std::tuple<MatrixXd,MatrixXd,MatrixXd,MatrixXd> get_lmc_function_weights(double x, double y, double z);
  std::tuple<std::vector<double>,std::vector<MatrixXd>> return_lmc_coefficients();
  void install_lmc_coefficients(std::vector<MatrixXd> tableau);

  std::tuple<MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd> get_disc_function_weights(double x, double y, double z);
  std::tuple<std::vector<double>,std::vector<MatrixXd>,std::vector<MatrixXd>> return_disc_coefficients();
  void install_disc_coefficients(std::vector<MatrixXd> costableau, std::vector<MatrixXd> sintableau);

};


MWLMC::MWLMC()
{
  initialise();
};

void MWLMC::initialise()
{
  // MW
  cout << "Initialising MW ... " << endl;
  MW = new SphExpansion(sph_cache_name_mw, model_file_mw, coef_file_mw, orient_file_mw);

  // LMC
  cout << "Initialising LMC ... " << endl;
  LMC = new SphExpansion(sph_cache_name_lmc, model_file_lmc, coef_file_lmc, orient_file_lmc);

  // MW
  cout << "Initialising MW disc ... " << endl;
  MWD = new CylExpansion(cyl_cache_name_mw, cyl_coef_name_mw, cyl_orient_name_mw);
};


void MWLMC::print_orbit(MatrixXd orbit, string orbitfile)
{
  ofstream outorbit;
  outorbit.open(orbitfile);

  outorbit << "# t [Gyr]; x [kpc]; y [kpc]; z [kpc]; vx [km/s] ; vy [km/s] ; vz [km/s] ; f_x [km/s/s] ; f_y [km/s/s] ; f_z [km/s/s];" << endl;

  for (int i=0; i<orbit.cols(); i++) {

    outorbit << setw(14) << orbit(9,i);

    for (int j=0; j<9; j++) {
      outorbit << setw(14) << orbit(j,i);
    }
    outorbit << endl;
  }

  outorbit.close();

}

std::vector<double> MWLMC::mwhalo_fields(double t, double x, double y, double z,
                                         //bool globalframe,
                                         int mwhharmonicflag, bool verbose)
{
  /*
    always comes out in the frame of the expansion: translate to inertial before input, if desired

   */

  // allocate workspace
   double tvir,xvir,yvir,zvir,xtmp,ytmp,ztmp;


  // translate time to virial
  physical_to_virial_time(t,tvir);

  /*
  vector<double> zerocoords(3),halozero(3);
   if (globalframe) {

     // get the zero of the disc: the coordinate centre
     return_centre(tvir, MWD->orient, zerocoords);
     // get the centre of the halo expansion
     return_centre(tvir, MW->orient, halozero);
     // we need to compute the requested position relative to the halo expansion.
     //   take the input, in the disc frame, and put in the halo frame

     // shift the expansion centres to the disc coordinate system
     xtmp = x - zerocoords[0];
     ytmp = y - zerocoords[1];
     ztmp = z - zerocoords[2];

   } else {
     xtmp = x;
     ytmp = y;
     ztmp = z;
   }
   */

   xtmp = x;
   ytmp = y;
   ztmp = z;

  // translate all times and positions into exp virial units
  physical_to_virial_length(xtmp,ytmp,ztmp,xvir,yvir,zvir);

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  tvir += reference_time;

  // get coefficients
  MatrixXd mwcoefs;
  MW->select_coefficient_time(tvir, mwcoefs);

  double rtmp,phitmp,thetatmp,r2tmp;
  double dens0,denstmp,tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys,pphys,dphys;

  // compute spherical coordinates in the frame of the MW expansion
  cartesian_to_spherical(xvir, yvir, zvir, rtmp, phitmp, thetatmp);

  // block from evaluating at tiny r
  if (rtmp < EEPS) rtmp = EEPS;

  // get all field values
  MW->determine_fields_at_point_sph(mwcoefs,
                                    rtmp,thetatmp,phitmp,
                                    dens0,denstmp,
                                    tpotl0,tpotl,
                                    fr,ft,fp,
                                    mwhharmonicflag);


  // for cases exactly along the z axis, blank out ft.
  if (verbose) {
    std::cout << std::setw(14) << rtmp
              << std::setw(14) << thetatmp
              << std::setw(14) << phitmp
              << std::setw(14) << xvir
              << std::setw(14) << yvir
              << std::setw(14) << zvir
              << std::setw(14) << fr
              << std::setw(14) << ft
              << std::setw(14) << fp
              << std::endl;
  }

  // for cases exactly along the z axis, blank out ft.
  if (isnan(ft)) {
    std::cout << "NaN ft" << std::endl;
    ft = 0.0;
  }

  // convert to cartesian
  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  if (verbose) {
    std::cout << "Cartesian: "
              << std::setw(14) << fxtmp
              << std::setw(14) << fytmp
              << std::setw(14) << fztmp
              << std::endl;
  }


  // translate all quantities to physical units
  virial_to_physical_density(denstmp,dphys);
  virial_to_physical_potential(tpotl,pphys);
  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // initialise output array
  std::vector<double> output(5);

  output[0] = fxphys;
  output[1] = fyphys;
  output[2] = fzphys;
  output[3] = dphys;
  output[4] = pphys;

  return output;

}



MatrixXd MWLMC::mwhalo_fields(double t, std::vector<double> x, std::vector<double> y, std::vector<double> z,
                              //bool globalframe,
                              int mwhharmonicflag, bool verbose)
{
  /*
    always comes out in the frame of the expansion: translate to inertial before input, if desired

   */

   double tvir;

   // translate all times and positions into exp virial units
   physical_to_virial_time(t,tvir);

   // reset time to have the correct system zero (e.g. pericentre is T=0)
   tvir += reference_time;

   // get coefficients
   MatrixXd mwcoefs;
   MW->select_coefficient_time(tvir, mwcoefs);

   double rtmp,phitmp,thetatmp,r2tmp;
   double dens0,denstmp,tpotl0,tpotl,fr,ft,fp;
   double fxtmp,fytmp,fztmp;

   double xphys,yphys,zphys,fxphys,fyphys,fzphys,pphys,dphys;

   double xvir,yvir,zvir,xtmp,ytmp,ztmp;

   // get the present-day MWD coordinates: the zero of the system
   vector<double> zerocoords(3);
   return_centre(reference_time, MWD->orient, zerocoords);

   // initialise output array
   int npositions = x.size();

   MatrixXd output;
   output.resize(npositions,5);

   for (int n=0;n<npositions;n++) {

/*
   if (globalframe) {

     // shift the expansion centres to the pericentre coordinate system
     xtmp = x[n] - zerocoords[0];
     ytmp = y[n] - zerocoords[1];
     ztmp = z[n] - zerocoords[2];

   } else {
     */
     xtmp = x[n];
     ytmp = y[n];
     ztmp = z[n];
   //}


  physical_to_virial_length(xtmp,ytmp,ztmp,xvir,yvir,zvir);

  // compute spherical coordinates in the frame of the MW expansion
  cartesian_to_spherical(xvir, yvir, zvir, rtmp, phitmp, thetatmp);

  // block from evaluating at tiny r
  if (rtmp < EEPS) rtmp = EEPS;

  // get all field values
  MW->determine_fields_at_point_sph(mwcoefs,
                                    rtmp,thetatmp,phitmp,
                                    dens0,denstmp,
                                    tpotl0,tpotl,
                                    fr,ft,fp,
                                    mwhharmonicflag);

  if (verbose) {
    std::cout << std::setw(14) << rtmp
              << std::setw(14) << thetatmp
              << std::setw(14) << phitmp
              << std::setw(14) << xvir
              << std::setw(14) << yvir
              << std::setw(14) << zvir
              << std::setw(14) << fr
              << std::setw(14) << ft
              << std::setw(14) << fp
              << std::endl;
  }

  // for cases exactly along the z axis, blank out ft.
  if (isnan(ft)) {
    if (verbose) std::cout << "NaN ft" << std::endl;
    ft = 0.0;
  }

  // convert to cartesian
  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);


                                if (verbose) {
                                  std::cout << "Cartesian: "
                                            << std::setw(14) << fxtmp
                                            << std::setw(14) << fytmp
                                            << std::setw(14) << fztmp
                                            << std::endl;
                                }

  // translate all quantities to physical units
  virial_to_physical_density(denstmp,dphys);
  virial_to_physical_potential(tpotl,pphys);
  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  output(n,0) = fxphys;
  output(n,1) = fyphys;
  output(n,2) = fzphys;
  output(n,3) = dphys;
  output(n,4) = pphys;
} // end npositions loop

  return output;

}




std::vector<double> MWLMC::lmc_fields(double t, double x, double y, double z,
                                      //bool globalframe,
                                      int lmcharmonicflag, bool verbose)
{
  /*
    always comes out in the frame of the expansion: translate to inertial before input, if desired

   */
   // translate all times and positions into exp virial units

  // translate all times and positions into exp virial units
  double tvir,xvir,yvir,zvir,xtmp,ytmp,ztmp;
  physical_to_virial_time(t,tvir);

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  tvir += reference_time;

  MatrixXd lmccoefs;
  LMC->select_coefficient_time(tvir, lmccoefs);

  double rtmp,phitmp,thetatmp,r2tmp;
  double dens0,denstmp,tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys,pphys,dphys;

/*
  if (globalframe) {

    // get the present-day MWD coordinates: the zero of the system
    vector<double> zerocoords(3);
    return_centre(reference_time, MWD->orient, zerocoords);

    // shift the expansion centres to the pericentre coordinate system
    xtmp = x - zerocoords[0];
    ytmp = y - zerocoords[1];
    ztmp = z - zerocoords[2];

  } else {*/
    xtmp = x;
    ytmp = y;
    ztmp = z;
  //}

  physical_to_virial_length(xtmp,ytmp,ztmp, xvir,yvir,zvir);

  // compute spherical coordinates in the frame of the LMC expansion
  cartesian_to_spherical(xvir, yvir, zvir, rtmp, phitmp, thetatmp);

  // block from evaluating at tiny r
  if (rtmp < EEPS) rtmp = EEPS;

  // get all field values
  LMC->determine_fields_at_point_sph(lmccoefs,
                                     rtmp,thetatmp,phitmp,
                                     dens0,denstmp,
                                     tpotl0,tpotl,
                                     fr,ft,fp,
                                     lmcharmonicflag);

   if (verbose) {
     std::cout << std::setw(14) << rtmp
               << std::setw(14) << thetatmp
               << std::setw(14) << phitmp
               << std::setw(14) << xvir
               << std::setw(14) << yvir
               << std::setw(14) << zvir
               << std::setw(14) << fr
               << std::setw(14) << ft
               << std::setw(14) << fp
               << std::endl;
   }

   // for cases exactly along the z axis, blank out ft.
   if (isnan(ft)) {
     if (verbose) std::cout << "NaN ft" << std::endl;
     ft = 0.0;
   }

  // convert to cartesian
  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  // translate all quantities to physical units
  virial_to_physical_density(denstmp,dphys);
  virial_to_physical_potential(tpotl,pphys);
  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // initialise output array
  std::vector<double> output(5);

  output[0] = fxphys;
  output[1] = fyphys;
  output[2] = fzphys;
  output[3] = dphys;
  output[4] = pphys;

  return output;

}


MatrixXd MWLMC::lmc_fields(double t, std::vector<double> x, std::vector<double> y, std::vector<double> z,
                           //bool globalframe,
                           int lmcharmonicflag, bool verbose)
{
  /*
    always comes out in the frame of the expansion: translate to inertial before input, if desired

   */
   // translate all times and positions into exp virial units
   vector<double> zerocoords(3);

  // translate all times and positions into exp virial units
  double tvir,xvir,yvir,zvir,xtmp,ytmp,ztmp;
  physical_to_virial_time(t,tvir);

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  tvir += reference_time;

  MatrixXd lmccoefs;
  LMC->select_coefficient_time(tvir, lmccoefs);

/*
  if (globalframe) {
    // get the present-day MWD coordinates: the zero of the system
    return_centre(reference_time, MWD->orient, zerocoords);
  }
  */

  double rtmp,phitmp,thetatmp,r2tmp;
  double dens0,denstmp,tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys,pphys,dphys;

  // initialise output array
  int npositions = x.size();

  MatrixXd output;
  output.resize(npositions,5);

  for (int n=0;n<npositions;n++) {
/*
    if (globalframe) {

      // shift the expansion centres to the pericentre coordinate system
      xtmp = x[n] - zerocoords[0];
      ytmp = y[n] - zerocoords[1];
      ztmp = z[n] - zerocoords[2];

    } else {*/
      xtmp = x[n];
      ytmp = y[n];
      ztmp = z[n];
  //  }

    physical_to_virial_length(xtmp,ytmp,ztmp, xvir,yvir,zvir);


    // compute spherical coordinates in the frame of the LMC expansion
    cartesian_to_spherical(xvir, yvir, zvir, rtmp, phitmp, thetatmp);

    // block from evaluating at tiny r
    if (rtmp < EEPS) rtmp = EEPS;

    // get all field values
    LMC->determine_fields_at_point_sph(lmccoefs,
                                       rtmp,thetatmp,phitmp,
                                       dens0,denstmp,
                                       tpotl0,tpotl,
                                       fr,ft,fp,
                                       lmcharmonicflag);


     if (verbose) {
       std::cout << std::setw(14) << rtmp
                 << std::setw(14) << thetatmp
                 << std::setw(14) << phitmp
                 << std::setw(14) << xvir
                 << std::setw(14) << yvir
                 << std::setw(14) << zvir
                 << std::setw(14) << fr
                 << std::setw(14) << ft
                 << std::setw(14) << fp
                 << std::endl;
     }

     // for cases exactly along the z axis, blank out ft.
     if (isnan(ft)) {
       if (verbose) std::cout << "NaN ft" << std::endl;
       ft = 0.0;
     }


    // convert to cartesian
    spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                  fr, fp, ft,
                                  fxtmp, fytmp, fztmp);

    // translate all quantities to physical units
    virial_to_physical_density(denstmp,dphys);
    virial_to_physical_potential(tpotl,pphys);
    virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

    output(n,0) = fxphys;
    output(n,1) = fyphys;
    output(n,2) = fzphys;
    output(n,3) = dphys;
    output(n,4) = pphys;
  } // end npositions loop


  return output;

}




std::vector<double> MWLMC::mwd_fields(double t, double x, double y, double z,
                                      //bool globalframe,
                                      int mwdharmonicflag,
                                      bool verbose)
{
  /*

   */

  // translate all times and positions into exp virial units
  double tvir,xvir,yvir,zvir;
  physical_to_virial_time(t,tvir);
  physical_to_virial_length(x,y,z, xvir,yvir,zvir);

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  tvir += reference_time;

  // get coefficients
  MatrixXd mwdcoscoefs,mwdsincoefs;
  MWD->select_coefficient_time(tvir, mwdcoscoefs, mwdsincoefs);
  //MWD->select_coefficient_time(0.0, mwdcoscoefs, mwdsincoefs);

  double phitmp,r2tmp;
  double tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  // the density fields are zombies right now: to be fixed with proper demand
  double dens0 = 0;
  double denstmp = 0;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys,pphys,dphys;

  // compute spherical coordinates in the frame of the MW expansion
  cartesian_to_cylindrical(xvir, yvir, r2tmp, phitmp);

  // same procedure for the disc
  MWD->determine_fields_at_point_cyl(mwdcoscoefs,mwdsincoefs,
                                     r2tmp,phitmp,zvir,
                                     tpotl0,tpotl,
                                     fr,fp,fztmp,mwdharmonicflag);

  cylindrical_forces_to_cartesian(r2tmp, phitmp,
                                  fr, fp,
                                  fxtmp, fytmp);

  // translate to physical units
  virial_to_physical_density(denstmp,dphys);
  virial_to_physical_potential(tpotl,pphys);
  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // return MW force

  // initialise output array
  std::vector<double> output(5);

  output[0] = fxphys;
  output[1] = fyphys;
  output[2] = fzphys;
  output[3] = dphys;
  output[4] = pphys;

  return output;

}



MatrixXd MWLMC::mwd_fields(double t, std::vector<double> x, std::vector<double> y, std::vector<double> z,
                                      //bool globalframe,
                                      int mwdharmonicflag,
                                      bool verbose)
{
  /*

   */

  // translate all times and positions into exp virial units
  double tvir,xvir,yvir,zvir;
  physical_to_virial_time(t,tvir);

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  tvir += reference_time;

  // get coefficients
  MatrixXd mwdcoscoefs,mwdsincoefs;
  MWD->select_coefficient_time(tvir, mwdcoscoefs, mwdsincoefs);

  double phitmp,r2tmp;
  double tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  // the density fields are zombies right now: to be fixed with proper demand
  double dens0 = 0;
  double denstmp = 0;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys,pphys,dphys;

  // initialise output array
  int npositions = x.size();

  MatrixXd output;
  output.resize(npositions,5);

  for (int n=0;n<npositions;n++) {

  physical_to_virial_length(x[n],y[n],z[n], xvir,yvir,zvir);

  // compute spherical coordinates in the frame of the MW expansion
  cartesian_to_cylindrical(xvir, yvir, r2tmp, phitmp);

  // same procedure for the disc
  MWD->determine_fields_at_point_cyl(mwdcoscoefs,mwdsincoefs,
                                     r2tmp,phitmp,zvir,
                                     tpotl0,tpotl,
                                     fr,fp,fztmp,mwdharmonicflag);

  cylindrical_forces_to_cartesian(r2tmp, phitmp,
                                  fr, fp,
                                  fxtmp, fytmp);

  // translate to physical units
  virial_to_physical_density(denstmp,dphys);
  virial_to_physical_potential(tpotl,pphys);
  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // return MW force

  output(n,0) = fxphys;
  output(n,1) = fyphys;
  output(n,2) = fzphys;
  output(n,3) = dphys;
  output(n,4) = pphys;
} // end npositions loop

  return output;

}




std::vector<double>  MWLMC::all_forces(double t, double x, double y, double z,
                                       //bool globalframe,
                                       int mwhharmonicflag, int mwdharmonicflag, int lmcharmonicflag,
                                       bool verbose)
{
  /*
  inputs are PHYSICAL
  t is Gyr
  x is kpc
  y is kpc
  z is kpc

  */
  // allocate virial workspace
  double tvir,xvir,yvir,zvir;

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  physical_to_virial_time(t,tvir);
  tvir += reference_time;

  // set up coordinate systems
  // initialise the centre vectors
  vector<double> zerocoords(3),mw_centre(3),lmc_centre(3),mwd_centre(3);

  // get the present-day MWD coordinates: the absolute zero of the system
  return_centre(reference_time, MWD->orient, zerocoords);

  //std::cout << "Coordinate zero:" << setw(14) << zerocoords[0] << setw(14) << zerocoords[1] << setw(14) << zerocoords[2] << endl;

  // get the centres of the expansions at the specified times in exp reference space
  return_centre(tvir,  MW->orient,  mw_centre);
  return_centre(tvir, LMC->orient, lmc_centre);
  return_centre(tvir, MWD->orient, mwd_centre);

  // zero out forces
  double fx = 0;
  double fy = 0;
  double fz = 0;

  // set up the output vector
  std::vector<double> output(3);
/*
  if (globalframe) {
    // we want to be in the GLOBAL frame of the MW disc: add the zero of the coordinate frame
    physical_to_virial_length(x + zerocoords[0],
                              y + zerocoords[1],
                              z + zerocoords[2],
                              xvir,yvir,zvir);
  } else {
    // we are in the inertial frame of the simulation: no adjustment needed
    */
    physical_to_virial_length(x,y,z,xvir,yvir,zvir);
  //}

  // get the initial coefficient values: the time here is in tvir units, so always start with 0
  MatrixXd tcoefsmw,tcoefslmc;
  MW->select_coefficient_time(tvir, tcoefsmw);
  LMC->select_coefficient_time(tvir, tcoefslmc);

  MatrixXd mwdcoscoefs,mwdsincoefs;
  MWD->select_coefficient_time(tvir, mwdcoscoefs, mwdsincoefs);

  // shift the expansion centres to the pericentre coordinate system
  for (int j=0;j<=2;j++) {
    mw_centre[j]  -= zerocoords[j];
    lmc_centre[j] -= zerocoords[j];
    mwd_centre[j] -= zerocoords[j];
  }

  if (verbose) {
    std::cout << "MWH centre (x,y,z)=(" << virial_to_physical_length(mw_centre[0])  << ","
                                        << virial_to_physical_length(mw_centre[1])  << ","
                                        << virial_to_physical_length(mw_centre[2])  << ")" << std::endl;
    std::cout << "MWD centre (x,y,z)=(" << virial_to_physical_length(mwd_centre[0]) << ","
                                        << virial_to_physical_length(mwd_centre[1]) << ","
                                        << virial_to_physical_length(mwd_centre[2]) << ")" << std::endl;
    std::cout << "LMC centre (x,y,z)=(" << virial_to_physical_length(lmc_centre[0]) << ","
                                        << virial_to_physical_length(lmc_centre[1]) << ","
                                        << virial_to_physical_length(lmc_centre[2]) << ")" << std::endl;
  }

  double rtmp,phitmp,thetatmp,r2tmp;
  double tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys;

  // compute spherical coordinates relative to the frame of the MW DISC expansion
  cartesian_to_spherical(xvir-mwd_centre[0], yvir-mwd_centre[1], zvir-mwd_centre[2], rtmp, phitmp, thetatmp);

  // block from evaluating at tiny r
  if (rtmp < EEPS) rtmp = EEPS;

  MW->determine_fields_at_point_sph(tcoefsmw,
                                    rtmp,thetatmp,phitmp,
                                    tpotl0,tpotl,
                                    fr,ft,fp,mwhharmonicflag);

  // for cases exactly along the z axis, blank out ft.
  if (isnan(fr)) {
    if (verbose) std::cout << "NaN fr MW" << std::endl;
    ft = 0.0;
  }   if (isnan(ft)) {
    if (verbose) std::cout << "NaN ft MW" << std::endl;
    ft = 0.0;
  }
  if (isnan(fp)) {
    if (verbose) std::cout << "NaN fp MW" << std::endl;
    fp = 0.0;
  }

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  r2tmp = sqrt((xvir-mwd_centre[0])*(xvir-mwd_centre[0]) + (yvir-mwd_centre[1])*(yvir-mwd_centre[1]));

  // block from evaluating at tiny r
  if (r2tmp < EEPS) r2tmp = EEPS;

  // same procedure for the disc
  MWD->determine_fields_at_point_cyl(mwdcoscoefs,mwdsincoefs,
                                     r2tmp,phitmp,zvir-mwd_centre[2],
                                     tpotl0,tpotl,
                                     fr,fp,fztmp,mwdharmonicflag);

  cylindrical_forces_to_cartesian(rtmp, phitmp,
                                  fr, fp,
                                  fxtmp, fytmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  // same procedure for LMC
  cartesian_to_spherical(xvir-lmc_centre[0], yvir-lmc_centre[1], zvir-lmc_centre[2], rtmp, phitmp, thetatmp);

  // block from evaluating at tiny r
  if (rtmp < EEPS) rtmp = EEPS;

  //cout << setw(14) << rtmp << setw(14) << phitmp << setw(14) << thetatmp << endl;

  LMC->determine_fields_at_point_sph(tcoefslmc,
                                     rtmp,thetatmp,phitmp,
                                     tpotl0,tpotl,
                                     fr,ft,fp,lmcharmonicflag);

   // for cases exactly along the z axis, blank out ft.
   if (isnan(fr)) {
     if (verbose) std::cout << "NaN fr LMC" << std::endl;
     ft = 0.0;
   }   if (isnan(ft)) {
     if (verbose) std::cout << "NaN ft LMC" << std::endl;
     ft = 0.0;
   }
   if (isnan(fp)) {
     if (verbose) std::cout << "NaN fp LMC" << std::endl;
     fp = 0.0;
   }

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  // reset to physical units
  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add LMC force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  output[0] = fx;
  output[1] = fy;
  output[2] = fz;

  return output;

}



MatrixXd  MWLMC::all_forces(double t, std::vector<double> x, std::vector<double> y, std::vector<double> z,
                            //bool globalframe,
                            int mwhharmonicflag, int mwdharmonicflag, int lmcharmonicflag,
                            bool verbose)
{
  /*

   */

  // zero out forces
  double fx = 0;
  double fy = 0;
  double fz = 0;

  // translate all times and positions into exp virial units
  double tvir,xvir,yvir,zvir;
  physical_to_virial_time(t,tvir);

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  tvir += reference_time;

  // get the initial coefficient values: the time here is in tvir units, so always start with 0
  MatrixXd tcoefsmw,tcoefslmc;
  MW->select_coefficient_time(tvir, tcoefsmw);
  LMC->select_coefficient_time(tvir, tcoefslmc);

  MatrixXd mwdcoscoefs,mwdsincoefs;
  MWD->select_coefficient_time(tvir, mwdcoscoefs, mwdsincoefs);

  // initialise the centre vectors
  vector<double> zerocoords(3),mw_centre(3),lmc_centre(3),mwd_centre(3);

  // get the present-day MWD coordinates: the zero of the system
  return_centre(reference_time, MWD->orient, zerocoords);

  //cout << "Coordinate zero:" << setw(14) << zerocoords[0] << setw(14) << zerocoords[1] << setw(14) << zerocoords[2] << endl;

  // get the centres of the expansions at the specified times in exp reference space
  return_centre(tvir,  MW->orient,  mw_centre);
  return_centre(tvir, LMC->orient, lmc_centre);
  return_centre(tvir, MWD->orient, mwd_centre);

  // shift the expansion centres to the pericentre coordinate system
  for (int j=0;j<=2;j++) {
    mw_centre[j]  -= zerocoords[j];
    lmc_centre[j] -= zerocoords[j];
    mwd_centre[j] -= zerocoords[j];
  }

  if (verbose) {
    cout << "MW virial centre (x,y,z)=(" << mw_centre[0] << ","<< mw_centre[1] << ","<< mw_centre[2] << ")" <<endl;
    cout << "MWD virial centre (x,y,z)=(" << mwd_centre[0] << ","<< mwd_centre[1] << ","<< mwd_centre[2] << ")" <<endl;
    cout << "LMC virial centre (x,y,z)=(" << lmc_centre[0] << ","<< lmc_centre[1] << ","<< lmc_centre[2] << ")" <<endl;
  }

  double rtmp,phitmp,thetatmp,r2tmp;
  double tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys;

  // initialise output array
  int npositions = x.size();

  MatrixXd output;
  output.resize(npositions,3);

  for (int n=0;n<npositions;n++) {

    physical_to_virial_length(x[n],y[n],z[n], xvir,yvir,zvir);

  // zero out forces
  fx = 0;
  fy = 0;
  fz = 0;

  // compute spherical coordinates in the frame of the MW expansion
  cartesian_to_spherical(xvir-mwd_centre[0], yvir-mwd_centre[1], zvir-mwd_centre[2], rtmp, phitmp, thetatmp);

  // block from evaluating at tiny r
  if (rtmp < EEPS) rtmp = EEPS;
  //cout << setw(14) << rtmp << setw(14) << phitmp << setw(14) << thetatmp << endl;

  MW->determine_fields_at_point_sph(tcoefsmw,
                                    rtmp,thetatmp,phitmp,
                                    tpotl0,tpotl,
                                    fr,ft,fp,mwhharmonicflag);

  // for cases exactly along the z axis, blank out ft.
  if (isnan(fr)) {
    if (verbose) std::cout << "NaN fr MW" << std::endl;
    ft = 0.0;
  }   if (isnan(ft)) {
    if (verbose) std::cout << "NaN ft MW" << std::endl;
    ft = 0.0;
  }
  if (isnan(fp)) {
    if (verbose) std::cout << "NaN fp MW" << std::endl;
    fp = 0.0;
  }

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  r2tmp = sqrt((xvir-mwd_centre[0])*(xvir-mwd_centre[0]) + (yvir-mwd_centre[1])*(yvir-mwd_centre[1]));

  // block from evaluating at tiny r
  if (r2tmp < EEPS) r2tmp = EEPS;

  // same procedure for the disc
  MWD->determine_fields_at_point_cyl(mwdcoscoefs,mwdsincoefs,
                                     r2tmp,phitmp,zvir-mwd_centre[2],
                                     tpotl0,tpotl,
                                     fr,fp,fztmp,mwdharmonicflag);

  cylindrical_forces_to_cartesian(rtmp, phitmp,
                                  fr, fp,
                                  fxtmp, fytmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  // same procedure for LMC
  cartesian_to_spherical(xvir-lmc_centre[0], yvir-lmc_centre[1], zvir-lmc_centre[2], rtmp, phitmp, thetatmp);

  // block from evaluating at tiny r
  if (rtmp < EEPS) rtmp = EEPS;

  //cout << setw(14) << rtmp << setw(14) << phitmp << setw(14) << thetatmp << endl;

  LMC->determine_fields_at_point_sph(tcoefslmc,
                                     rtmp,thetatmp,phitmp,
                                     tpotl0,tpotl,
                                     fr,ft,fp,lmcharmonicflag);

   // for cases exactly along the z axis, blank out ft.
   if (isnan(fr)) {
     if (verbose) std::cout << "NaN fr LMC" << std::endl;
     ft = 0.0;
   }   if (isnan(ft)) {
     if (verbose) std::cout << "NaN ft LMC" << std::endl;
     ft = 0.0;
   }
   if (isnan(fp)) {
     if (verbose) std::cout << "NaN fp LMC" << std::endl;
     fp = 0.0;
   }

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  // reset to physical units
  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add LMC force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  output(n,0) = fx;
  output(n,1) = fy;
  output(n,2) = fz;
}

  return output;

}


MatrixXd MWLMC::get_trajectories(double dt, bool virial)
{
  // get the trajectories with the specified time sampling (returns native time sampling if not specified)
  int nint = reference_time/dt;

  MatrixXd trajectories;
  trajectories.resize(nint,13);

  double tphys,xphys;


  for (int n=0;n<nint;n++)
  {
    // always start in native time, as we know the simulation bounds
    std::vector<double> trajectorytmp = get_expansion_centres_virial(n*dt);

    if (virial) {

      trajectories(n,0) = n*dt;
      for (int j=0;j<12;j++) trajectories(n,j+1) = trajectorytmp[j];

    } else {

      // translate back to physical units
      virial_to_physical_time(n*dt,tphys);
      trajectories(n,0) = tphys;

      for (int j=0;j<12;j++) {
        trajectories(n,j+1) = virial_to_physical_length(trajectorytmp[j]);
      }

    }
  }

  return trajectories;

}

MatrixXd MWLMC::get_lmc_trajectory(double rewindtime, double dt)
{
  // get the trajectories with the specified time sampling (returns native time sampling if not specified)

  double dtvir,rewindtimevir,tphys,xphys;

  // convert times to virial for under the hood operation
  physical_to_virial_time(dt,dtvir);
  physical_to_virial_time(rewindtime,rewindtimevir);

  // calculate the number of steps to take
  int nint = rewindtimevir/dtvir;

  MatrixXd trajectory;
  trajectory.resize(nint,4);


  for (int n=0;n<nint;n++)
  {
    // always start in native time, as we know the simulation bounds
    std::vector<double> trajectorytmp = get_lmc_centre_virial(reference_time-n*dtvir,false);

    // translate time back to physical units
    virial_to_physical_time(-(n*dtvir),tphys);

    // assign time
    trajectory(n,0) = tphys;

    // translate lengths back to physical units
    for (int j=0;j<3;j++) {
      trajectory(n,j+1) = virial_to_physical_length(trajectorytmp[j]);
    }
  }

  return trajectory;

}


std::vector<double> MWLMC::get_expansion_centres_physical(double t, bool verbose)
{
  // check for a valid time
  if (t > 0.0) {
    std::cout << "Cannot select a time after the present day! Setting to present day..." << std::endl;
    t = 0.0;
  }

  // make one big return container
  std::vector<double> centres(12);

  // translate all times and positions into exp virial units
  double tvir;
  physical_to_virial_time(t,tvir);

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  tvir += reference_time;

  centres = get_expansion_centres_virial(tvir, verbose);

  // recast to physical distances
  for (int i=0;i<12;i++) centres[i] = virial_to_physical_length(centres[i]);

  return centres;

}

void MWLMC::reset_mw_coefficients()
{

  MW->reset_coefficients();

}

void MWLMC::reset_all_coefficients()
{

  MW->reset_coefficients();
  LMC->reset_coefficients();
  MWD->reset_coefficients();

}

std::tuple<MatrixXd,MatrixXd,MatrixXd,MatrixXd> MWLMC::get_mw_function_weights(double x, double y, double z)
{

  // call translator to virial units
  double xvir,yvir,zvir;
  physical_to_virial_length(x,y,z, xvir,yvir,zvir);

  // translate to spherical coordinates
  double rtmp,phitmp,thetatmp;
  cartesian_to_spherical(xvir, yvir, zvir, rtmp, phitmp, thetatmp);

  std::tuple<MatrixXd,MatrixXd,MatrixXd,MatrixXd> X;
  X = MW->determine_weights_at_point_sph(rtmp,thetatmp,phitmp);

  MatrixXd pot,fr,ft,fp;
  pot  = std::get<0>(X);
  fr   = std::get<1>(X);
  ft   = std::get<2>(X);
  fp   = std::get<3>(X);

  // convert forces
  MatrixXd fxtmp,fytmp,fztmp;
  spherical_forces_to_cartesian(rtmp,phitmp,thetatmp,fr,fp,ft,fxtmp,fytmp,fztmp);

  // rescale potential and forces to physical units
  MatrixXd potphys,fxphys,fyphys,fzphys;
  virial_to_physical_force(fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);
  virial_to_physical_potential(pot, potphys);

  return make_tuple(potphys,fxphys,fyphys,fzphys);
}

std::tuple<std::vector<double>,std::vector<MatrixXd>> MWLMC::return_mw_coefficients()
{

  std::tuple<std::vector<double>,std::vector<MatrixXd>> X;
  X = MW->return_coefficients();
  std::vector<double>   virialtimes = std::get<0>(X);
  std::vector<MatrixXd> coeftableau = std::get<1>(X);

  // convert the virial times to physical times
  std::vector<double> physicaltimes;
  physicaltimes.resize(virialtimes.size());
  for (int i=0;i<virialtimes.size();i++) {
    physicaltimes[i] = virial_to_physical_time_return(virialtimes[i] - reference_time);
  }

  return make_tuple(physicaltimes,coeftableau);
}

void MWLMC::install_mw_coefficients(std::vector<MatrixXd> tableau)
{
  MW->install_coefficients(tableau);
}

std::tuple<MatrixXd,MatrixXd,MatrixXd,MatrixXd> MWLMC::get_lmc_function_weights(double x, double y, double z)
{

  // call translator to virial units
  double xvir,yvir,zvir;
  physical_to_virial_length(x,y,z, xvir,yvir,zvir);

  // translate to spherical coordinates
  double rtmp,phitmp,thetatmp;
  cartesian_to_spherical(xvir, yvir, zvir, rtmp, phitmp, thetatmp);

  std::tuple<MatrixXd,MatrixXd,MatrixXd,MatrixXd> X;
  X = LMC->determine_weights_at_point_sph(rtmp,thetatmp,phitmp);

  MatrixXd pot,fr,ft,fp;
  pot  = std::get<0>(X);
  fr   = std::get<1>(X);
  ft   = std::get<2>(X);
  fp   = std::get<3>(X);

  // convert forces
  MatrixXd fxtmp,fytmp,fztmp;
  spherical_forces_to_cartesian(rtmp,phitmp,thetatmp,fr,fp,ft,fxtmp,fytmp,fztmp);

  // rescale potential and forces to physical units
  MatrixXd potphys,fxphys,fyphys,fzphys;
  virial_to_physical_force(fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);
  virial_to_physical_potential(pot, potphys);

  return make_tuple(potphys,fxphys,fyphys,fzphys);

}

std::tuple<std::vector<double>,std::vector<MatrixXd>>  MWLMC::return_lmc_coefficients()
{
  std::tuple<std::vector<double>,std::vector<MatrixXd>> X;
  X = LMC->return_coefficients();
  std::vector<double> virialtimes   = std::get<0>(X);
  std::vector<MatrixXd> coeftableau = std::get<1>(X);

  // convert the virial times to physical times
  std::vector<double> physicaltimes;
  physicaltimes.resize(virialtimes.size());
  for (int i=0;i<virialtimes.size();i++) {
    physicaltimes[i] = virial_to_physical_time_return(virialtimes[i] - reference_time);
  }

  return make_tuple(physicaltimes,coeftableau);

}

void MWLMC::install_lmc_coefficients(std::vector<MatrixXd> tableau)
{
  LMC->install_coefficients(tableau);
}

std::tuple<MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd> MWLMC::get_disc_function_weights(double x, double y, double z)
{

  // call translator to virial units
  double xvir,yvir,zvir;
  physical_to_virial_length(x,y,z, xvir,yvir,zvir);

  // compute spherical coordinates in the frame of the expansion
  double r2tmp,phitmp;
  cartesian_to_cylindrical(xvir,yvir,r2tmp,phitmp);
  // offset to the correct expansion centre (if these are coming in inertial space? or should the user have to do that?)

  // call out for weights
  std::tuple<MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd,MatrixXd> X;
  X = MWD->determine_weights_at_point_cyl(r2tmp,phitmp,zvir);

  // unpack
  MatrixXd potc,frc,fpc,fzc,pots,frs,fps,fzs;
  potc = std::get<0>(X);
  frc  = std::get<1>(X);
  fpc  = std::get<2>(X);
  fzc  = std::get<3>(X);
  pots = std::get<4>(X);
  frs  = std::get<5>(X);
  fps  = std::get<6>(X);
  fzs  = std::get<7>(X);

  // convert forces
  MatrixXd fxtmpc,fytmpc,fxtmps,fytmps;
  cylindrical_forces_to_cartesian(r2tmp, phitmp,
                                  frc, fpc,
                                  fxtmpc, fytmpc);

  cylindrical_forces_to_cartesian(r2tmp, phitmp,
                                  frs, fps,
                                  fxtmps, fytmps);

  // rescale potential and forces to physical units
  MatrixXd potcphys,fxphysc,fyphysc,fzphysc,potsphys,fxphyss,fyphyss,fzphyss;
  virial_to_physical_force(fxtmpc,fytmpc,fzc,fxphysc,fyphysc,fzphysc);
  virial_to_physical_force(fxtmps,fytmps,fzs,fxphyss,fyphyss,fzphyss);
  virial_to_physical_potential(potc, potcphys);
  virial_to_physical_potential(pots, potsphys);

  // return
  return make_tuple(potcphys,fxphysc,fyphysc,fzphysc,potsphys,fxphyss,fyphyss,fzphyss);

}

std::tuple<std::vector<double>,std::vector<MatrixXd>,std::vector<MatrixXd>>  MWLMC::return_disc_coefficients()
{
  std::tuple<std::vector<double>,std::vector<MatrixXd>,std::vector<MatrixXd>> X;
  X = MWD->return_coefficients();
  std::vector<double> virialtimes  = std::get<0>(X);
  std::vector<MatrixXd> costableau = std::get<1>(X);
  std::vector<MatrixXd> sintableau = std::get<2>(X);

  // convert the virial times to physical times
  std::vector<double> physicaltimes;
  physicaltimes.resize(virialtimes.size());
  for (int i=0;i<virialtimes.size();i++) {
    physicaltimes[i] = virial_to_physical_time_return(virialtimes[i] - reference_time);
  }

  return make_tuple(physicaltimes,costableau,sintableau);
}

void MWLMC::install_disc_coefficients(std::vector<MatrixXd> costableau, std::vector<MatrixXd> sintableau)
{
  MWD->install_coefficients(costableau, sintableau);
}


std::vector<double> MWLMC::get_expansion_centres_virial(double tvir, bool verbose)
{
  // check for a valid time
  if (tvir > reference_time) {
    std::cout << "Cannot select a time after the present day! Setting to present day..." << std::endl;
    tvir = reference_time;
  }

  // make one big return container
  std::vector<double> centres(12);

  // initialise temporary centre vectors
  vector<double> zerocoords(3),mw_centre(3),lmc_centre(3),mwd_centre(3);

  // get the present-day MWD coordinates: the zero of the total
  return_centre(tvir, MWD->orient, zerocoords);

  if (verbose) std::cout << "Coordinate zero:" << setw(14) << zerocoords[0] << setw(14) << zerocoords[1] << setw(14) << zerocoords[2] << std::endl;

  // get the centres of the expansions at the specified times in exp reference space
  return_centre(tvir,  MW->orient,  mw_centre);
  return_centre(tvir, LMC->orient, lmc_centre);
  return_centre(tvir, MWD->orient, mwd_centre);

  int indx;
  double dt;
  find_time_index(tvir, LMC->orient, indx, dt);

  // fill return vector
  for (int i=0;i<3;i++) centres[i  ] = zerocoords[i];
  for (int i=0;i<3;i++) centres[i+3] = mw_centre[i];
  for (int i=0;i<3;i++) centres[i+6] = lmc_centre[i];
  for (int i=0;i<3;i++) centres[i+9] = mwd_centre[i];

  return centres;

}



std::vector<double> MWLMC::get_expansion_centre_velocities_physical(double t, bool verbose)
{
  // check for a valid time
  if (t > 0.0) {
    std::cout << "Cannot select a time after the present day! Setting to present day..." << std::endl;
    t = 0.0;
  }

  // make one big return container
  std::vector<double> vcentres(12);

  // translate all times and positions into exp virial units
  double tvir;
  physical_to_virial_time(t,tvir);

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  tvir += reference_time;

  vcentres = get_expansion_centre_velocities_virial(tvir, verbose);

  // recast to physical velocities
  for (int i=0;i<12;i++) vcentres[i] = virial_to_physical_velocity(vcentres[i]);

  return vcentres;

}



std::vector<double> MWLMC::get_expansion_centre_velocities_virial(double tvir, bool verbose)
{
  // check for a valid time
  if (tvir > reference_time) {
    std::cout << "Cannot select a time after the present day! Setting to present day..." << std::endl;
    tvir = reference_time;
  }

  // make one big return container
  std::vector<double> vcentres(12);

  // initialise temporary centre vectors
  vector<double> zerovcoords(3),mw_vcentre(3),lmc_vcentre(3),mwd_vcentre(3);

  // get the present-day MWD coordinates: the zero of the total
  return_vel_centre(tvir, MWD->orient, zerovcoords);

  if (verbose) std::cout << "Coordinate zero:" << setw(14) << zerovcoords[0] << setw(14) << zerovcoords[1] << setw(14) << zerovcoords[2] << std::endl;

  // get the centres of the expansions at the specified times in exp reference space
  return_vel_centre(tvir,  MW->orient,  mw_vcentre);
  return_vel_centre(tvir, LMC->orient, lmc_vcentre);
  return_vel_centre(tvir, MWD->orient, mwd_vcentre);

  int indx;
  double dt;
  find_time_index(tvir, LMC->orient, indx, dt);

  // fill return vector
  for (int i=0;i<3;i++) vcentres[i  ] = zerovcoords[i];
  for (int i=0;i<3;i++) vcentres[i+3] = mw_vcentre[i];
  for (int i=0;i<3;i++) vcentres[i+6] = lmc_vcentre[i];
  for (int i=0;i<3;i++) vcentres[i+9] = mwd_vcentre[i];

  return vcentres;

}


std::vector<double> MWLMC::get_lmc_centre_virial(double tvir, bool verbose)
{
  /*
  get_lmc_centre_virial
    computes the location of the LMC centre relative to the MW disc

  */
  // check for a valid time
  if (tvir > reference_time) {
    std::cout << "Cannot select a time after the present day! Setting to present day..." << std::endl;
    tvir = reference_time;
  }

  // make one big return container
  std::vector<double> centres(3);

  // initialise temporary centre vectors
  vector<double> zerocoords(3),lmc_centre(3);

  // get the present-day MWD coordinates: the zero of the total
  return_centre(tvir, MWD->orient, zerocoords);

  // get the centres of the expansions at the specified times in exp reference space
  return_centre(tvir, LMC->orient, lmc_centre);

  // fill return vector
  for (int i=0;i<3;i++) centres[i  ] = lmc_centre[i] - zerocoords[i];

  return centres;

}


void MWLMC::all_forces_coefs(MatrixXd mwcoefs, MatrixXd lmccoefs,
                             MatrixXd mwdcoscoefs, MatrixXd mwdsincoefs,
                             double t, double x, double y, double z,
                             double& fx, double& fy, double& fz,
                             int mwhharmonicflag, int mwdharmonicflag, int lmcharmonicflag,
                             bool verbose)
{
  /*
    specs: take a time, x,y,z; return x,y,z forces, in physical units
    input/output units must be physical
    where (x,y,z)=(0,0,0) is the present-day DISC galactic centre
       and      t=0 is the present day (previous times are negative)

   */

  // zero out any residual force values
  fx = fy = fz = 0.;

  // translate all times and positions into exp virial units
  double tvir,xvir,yvir,zvir;
  physical_to_virial_time(t,tvir);
  physical_to_virial_length(x,y,z, xvir,yvir,zvir);

  // reset time to have the correct system zero (e.g. pericentre is T=0)
  tvir += reference_time;

  // initialise the centre vectors
  vector<double> zerocoords(3),mw_centre(3),lmc_centre(3),mwd_centre(3);
  return_centre(reference_time, MWD->orient, zerocoords); // the absolute zero of the coordinate system
  return_centre(tvir,  MW->orient,  mw_centre);
  return_centre(tvir, LMC->orient, lmc_centre);
  return_centre(tvir, MWD->orient, mwd_centre);

  // shift the expansion centres to the pericentre coordinate system
  for (int j=0;j<=2;j++) {
    mw_centre[j]  -= zerocoords[j];
    lmc_centre[j] -= zerocoords[j];
    mwd_centre[j] -= zerocoords[j];
  }

  if (verbose) {
    cout << "MW virial centre (x,y,z)=(" << mw_centre[0] << ","<< mw_centre[1] << ","<< mw_centre[2] << ")" <<endl;
    cout << "MWD virial centre (x,y,z)=(" << mwd_centre[0] << ","<< mwd_centre[1] << ","<< mwd_centre[2] << ")" <<endl;
    cout << "LMC virial centre (x,y,z)=(" << lmc_centre[0] << ","<< lmc_centre[1] << ","<< lmc_centre[2] << ")" <<endl;
  }

  double rtmp,phitmp,thetatmp,r2tmp;
  double tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys;

  // compute spherical coordinates in the frame of the MW expansion
  cartesian_to_spherical(xvir-mw_centre[0], yvir-mw_centre[1], zvir-mw_centre[2], rtmp, phitmp, thetatmp);

  //cout << setw(14) << rtmp << setw(14) << phitmp << setw(14) << thetatmp << endl;

  MW->determine_fields_at_point_sph(mwcoefs,
                                    rtmp,thetatmp,phitmp,
                                    tpotl0,tpotl,
                                    fr,ft,fp,mwhharmonicflag);

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  r2tmp = sqrt((xvir-mwd_centre[0])*(xvir-mwd_centre[0]) + (yvir-mwd_centre[1])*(yvir-mwd_centre[1]));

  // same procedure for the disc
  MWD->determine_fields_at_point_cyl(mwdcoscoefs,mwdsincoefs,
                                     r2tmp,phitmp,zvir-mwd_centre[2],
                                     tpotl0,tpotl,
                                     fr,fp,fztmp,mwdharmonicflag);

  cylindrical_forces_to_cartesian(rtmp, phitmp,
                                  fr, fp,
                                  fxtmp, fytmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  // same procedure for LMC
  cartesian_to_spherical(xvir-lmc_centre[0], yvir-lmc_centre[1], zvir-lmc_centre[2], rtmp, phitmp, thetatmp);

  //cout << setw(14) << rtmp << setw(14) << phitmp << setw(14) << thetatmp << endl;

  LMC->determine_fields_at_point_sph(lmccoefs,
                                     rtmp,thetatmp,phitmp,
                                     tpotl0,tpotl,
                                     fr,ft,fp,lmcharmonicflag);

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  // reset to physical units
  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add LMC force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

}


void MWLMC::mw_forces_coefs(MatrixXd mwcoefs,
                            MatrixXd mwdcoscoefs, MatrixXd mwdsincoefs,
                            double t, double x, double y, double z,
                            double& fx, double& fy, double& fz,
                            int mwhharmonicflag, int mwdharmonicflag,
                            bool verbose)
{
  /*
    No translation: will only ever happen in the local frame of the MW.

   @IMPROVE: t is a dead parameter, can remove from here and elsewhere
   */

  // zero out forces
  fx = fy = fz = 0.;

  // translate all times and positions into exp virial units
  double xvir,yvir,zvir;
  physical_to_virial_length(x,y,z, xvir,yvir,zvir);

  double rtmp,phitmp,thetatmp,r2tmp;
  double tpotl0,tpotl,fr,ft,fp;
  double fxtmp,fytmp,fztmp;

  double xphys,yphys,zphys,fxphys,fyphys,fzphys;

  // compute spherical coordinates
  cartesian_to_spherical(xvir, yvir, zvir, rtmp, phitmp, thetatmp);

  MW->determine_fields_at_point_sph(mwcoefs,
                                    rtmp,thetatmp,phitmp,
                                    tpotl0,tpotl,
                                    fr,ft,fp,mwhharmonicflag);

  spherical_forces_to_cartesian(rtmp, phitmp, thetatmp,
                                fr, fp, ft,
                                fxtmp, fytmp, fztmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

  r2tmp = sqrt(xvir*xvir + yvir*yvir);

  // same procedure for the disc
  MWD->determine_fields_at_point_cyl(mwdcoscoefs,mwdsincoefs,
                                     r2tmp,phitmp,zvir,
                                     tpotl0,tpotl,
                                     fr,fp,fztmp,mwdharmonicflag);

  cylindrical_forces_to_cartesian(rtmp, phitmp,
                                  fr, fp,
                                  fxtmp, fytmp);

  virial_to_physical_force (fxtmp,fytmp,fztmp,fxphys,fyphys,fzphys);

  // add MW force to total
  fx += fxphys;
  fy += fyphys;
  fz += fzphys;

}



MatrixXd MWLMC::mworbit(vector<double> xinit,
                        vector<double> vinit,
                        double tbegin,
                        double tend,
                        double dt)
{
  /*
  integration of just the Milky Way components
   */

  // allocate workspace
  double fx,fy,fz,tvir;

  // allocate output array
  MatrixXd orbit;

  // number of integration steps to take
  int nint = (int)((tend - tbegin)/dt);


  if (nint<0) {
    std::cerr << "orbit: ending time must be before beginning time. Exiting." << std::endl;
    exit(-1);
  }

  // include the forces for now
  orbit.resize(10,nint);

   // set times
   double tvirbegin  = physical_to_virial_time_return(tbegin) + reference_time;
   double dtvir      = physical_to_virial_time_return(dt);
   double tphysbegin = tbegin;
   double dtphys     = dt;

   orbit(0,0) = xinit[0];
   orbit(1,0) = xinit[1];
   orbit(2,0) = xinit[2];
   orbit(3,0) = vinit[0];
   orbit(4,0) = vinit[1];
   orbit(5,0) = vinit[2];
   // 6-8 are forces, set below...
   orbit(9,0) = tphysbegin;

   //now step forward one, using leapfrog (drift-kick-drift) integrator?
   //    https://en.wikipedia.org/wiki/Leapfrog_integration
   //
   int step = 1;


   if (tvirbegin < 0) {
     std::cout << "fullintegrate.h: tvirbegin is less than 0! --> " << tvirbegin << std::endl;
   }

   // get the initial coefficient values: the time here is in tvir units, so always start with 0
   MatrixXd tcoefsmw,tcoefslmc, mwcoscoefs,mwsincoefs;
   MW->select_coefficient_time(0., tcoefsmw);
   MWD->select_coefficient_time(0., mwcoscoefs, mwsincoefs);

   // return forces for the initial step
   mw_forces_coefs(tcoefsmw, mwcoscoefs, mwsincoefs,
                   0.0, orbit(0,0),orbit(1,0),orbit(2,0),
                   fx, fy, fz,
                   127,127);

   orbit(6,0) = fx;
   orbit(7,0) = fy;
   orbit(8,0) = fz;

   int j;
   double tvirnow, tphysnow;

   for (step=1; step<nint; step++) {

     // advance timestep: this is in virial units by definition.
     tvirnow  = tvirbegin  + dtvir*step;
     tphysnow = tphysbegin + dtphys*step;
     orbit(9,step) = tphysnow; // record the time

     // advance positions
     for (j=0; j<3; j++) {
       orbit(j,step) = orbit(j,step-1)   + (orbit(j+3,step-1)*dt  )  + (0.5*orbit(j+6,step-1)  * (dt*dt));
     }

     // this call goes out in physical units
     mw_forces_coefs(tcoefsmw, mwcoscoefs, mwsincoefs,
                     0.0, orbit(0,step),orbit(1,step),orbit(2,step),
                     fx, fy, fz,
                     127,127);

     orbit(6,step) = fx;
     orbit(7,step) = fy;
     orbit(8,step) = fz;

     // advance velocities
     for (j=3; j<6; j++) {
       orbit(j,step) = orbit(j,step-1) + (0.5*(orbit(j+3,step-1)+orbit(j+3,step))  * dt );
     }

   }

   return orbit;

 }



std::vector< MatrixXd > MWLMC::mworbit(MatrixXd xinit,
                                       MatrixXd vinit,
                                       double tbegin,
                                       double tend,
                                       double dt)
 {
   /*
   integration of just the Milky Way components
    */

    double fx,fy,fz,tvir;

    int norbits = xinit.rows();
    std::vector< MatrixXd > orbit;
    orbit.resize(norbits);

    // number of integration steps to take
    int nint = (int)((tend - tbegin)/dt);


    if (nint<0) {
      std::cerr << "orbit: ending time must be before beginning time. Exiting." << std::endl;
      exit(-1);
    }

    // set times
    double tvirbegin  = physical_to_virial_time_return(tbegin) + reference_time;
    double dtvir      = physical_to_virial_time_return(dt);
    double tphysbegin = tbegin;
    double dtphys     = dt;

    // initialise values
    for (int n=0;n<norbits;n++) {
      orbit[n].resize(10,nint);
      orbit[n](0,0) = xinit(n,0);
      orbit[n](1,0) = xinit(n,1);
      orbit[n](2,0) = xinit(n,2);
      orbit[n](3,0) = vinit(n,0);
      orbit[n](4,0) = vinit(n,1);
      orbit[n](5,0) = vinit(n,2);
      // 6-8 are forces, set below...
      orbit[n](9,0) = tphysbegin;
    }


    //now step forward one, using leapfrog (drift-kick-drift) integrator?
    //    https://en.wikipedia.org/wiki/Leapfrog_integration
    //
    int step = 1;


    if (tvirbegin < 0) {
      std::cout << "fullintegrate.h: tvirbegin is less than 0! --> " << tvirbegin << std::endl;
    }

    // get the initial coefficient values: the time here is in tvir units, so always start with 0
    MatrixXd tcoefsmw,tcoefslmc, mwcoscoefs,mwsincoefs;
    MW->select_coefficient_time(0., tcoefsmw);
    MWD->select_coefficient_time(0., mwcoscoefs, mwsincoefs);

    for (int n=0;n<norbits;n++) {
      // return forces for the initial step
      mw_forces_coefs(tcoefsmw, mwcoscoefs, mwsincoefs,
                      0., orbit[n](0,0),orbit[n](1,0),orbit[n](2,0),
                      fx, fy, fz,
                      127,127);

      orbit[n](6,0) = fx;
      orbit[n](7,0) = fy;
      orbit[n](8,0) = fz;

    }

    int j;
    double tvirnow, tphysnow;

    for (step=1; step<nint; step++) {

      // advance timestep: this is in virial units by definition.
      tvirnow  = tvirbegin  + dtvir*step;
      tphysnow = tphysbegin + dtphys*step;


      for (int n=0;n<norbits;n++) {

        // record the time
        orbit[n](9,step) = tphysnow;

        // advance positions
        for (j=0; j<3; j++) {
          orbit[n](j,step) = orbit[n](j,step-1)   + (orbit[n](j+3,step-1)*dt  )  + (0.5*orbit[n](j+6,step-1)  * (dt*dt));
        }

        // this call goes out in physical units
        mw_forces_coefs(tcoefsmw, mwcoscoefs, mwsincoefs,
                        0., orbit[n](0,step),orbit[n](1,step),orbit[n](2,step),
                        fx, fy, fz,
                        127,127);

        orbit[n](6,step) = fx;
        orbit[n](7,step) = fy;
        orbit[n](8,step) = fz;

        // advance velocities
        for (j=3; j<6; j++) {
          orbit[n](j,step) = orbit[n](j,step-1) + (0.5*(orbit[n](j+3,step-1)+orbit[n](j+3,step))  * dt );
        }
      } // norbit loop

    } // step loop

    return orbit;

}



MatrixXd MWLMC::rewind(vector<double> xinit,
                       vector<double> vinit,
                       double dt,
                       int mwhharmonicflag, int mwdharmonicflag, int lmcharmonicflag,
                       double rewind_time,
                       bool discframe,
                       bool verbose)
{
 /*
 rewind
 -------------
 for a given position and velocity, rewind the orbit by a specified amount of time

 takes physical units

  */

  // step 0: allocate workspace
  double fx,fy,fz,tvir;

  // step 1: get the coordinate frame of the MW disc: the absolute frame of the system
  vector<double> disccoords(3),discvelcoords(3);
  vector<double> initcoords(3),initvelcoords(3);
  return_centre(reference_time, MWD->orient, initcoords);
  return_vel_centre(reference_time, MWD->orient, initvelcoords);

  if (verbose) {
    std::cout << setw(14) << reference_time
              << setw(14) << virial_to_physical_length(initcoords[0])
              << setw(14) << virial_to_physical_length(initcoords[1])
              << setw(14) << virial_to_physical_length(initcoords[2]) << std::endl;
    std::cout << setw(14) << "(u,v,w)"
              << setw(14) << virial_to_physical_velocity(initvelcoords[0])
              << setw(14) << virial_to_physical_velocity(initvelcoords[1])
              << setw(14) << virial_to_physical_velocity(initvelcoords[2]) << std::endl;
  }

  // step 2: set beginning times and timesteps
  double tvirbegin  = reference_time;
  double dtvir      = physical_to_virial_time_return(dt);
  double tphysbegin = 0.;
  double dtphys     = dt;

  // number of integration steps to take
  int nint = (int)(rewind_time/dt);

  if (nint<0) {
    std::cerr << "rewind: rewind_time must be positive. Exiting." << std::endl;
    exit(-1);
  }

  // make the output frame
  MatrixXd orbit;
  orbit.resize(10,nint);

  // initialise positions and (inverted) velocities for backwards integration
  orbit(0,0) =  xinit[0];
  orbit(1,0) =  xinit[1];
  orbit(2,0) =  xinit[2];
  orbit(3,0) = -vinit[0];
  orbit(4,0) = -vinit[1];
  orbit(5,0) = -vinit[2];
  //6-8 are forces, set below...
  orbit(9,0) = tphysbegin;

  //now step forward one, using leapfrog (drift-kick-drift) integrator
  //    https://en.wikipedia.org/wiki/Leapfrog_integration
  //
  int step = 1;

  // get the initial coefficient values: the time here is in tvir units, so always start with 0
  MatrixXd tcoefsmw,tcoefslmc, mwcoscoefs,mwsincoefs;
  MW->select_coefficient_time(tvirbegin, tcoefsmw);
  LMC->select_coefficient_time(tvirbegin, tcoefslmc);
  MWD->select_coefficient_time(tvirbegin, mwcoscoefs, mwsincoefs);

  // return forces for the initial step.
  all_forces_coefs(tcoefsmw, tcoefslmc, mwcoscoefs, mwsincoefs,
                  tphysbegin, orbit(0,0),orbit(1,0),orbit(2,0),
                  fx, fy, fz,
                  mwhharmonicflag, mwdharmonicflag, lmcharmonicflag);

  orbit(6,0) = fx;
  orbit(7,0) = fy;
  orbit(8,0) = fz;

  int j;
  double tvirnow, tphysnow;

  for (step=1; step<nint; step++) {

    // 'advance' timestep: this is in virial units by definition.
    tvirnow  = tvirbegin  - dtvir*step;
    tphysnow = tphysbegin - dtphys*step;
    orbit(9,step) = tphysnow; // record the time

    MW->select_coefficient_time(tvirnow, tcoefsmw);
    LMC->select_coefficient_time(tvirnow, tcoefslmc);
    MWD->select_coefficient_time(tvirnow, mwcoscoefs, mwsincoefs);

    // 'advance' positions
    for (j=0; j<3; j++) {
      orbit(j,step) = orbit(j,step-1)   + (orbit(j+3,step-1)*dt  )  + (0.5*orbit(j+6,step-1)  * (dt*dt));
    }

    // this call goes out in PHYSICAL units
    all_forces_coefs(tcoefsmw, tcoefslmc, mwcoscoefs, mwsincoefs,
                     // need to shift the time by one unit to get the proper centre
                     tphysnow - dtphys, orbit(0,step),orbit(1,step),orbit(2,step),
                     fx, fy, fz,
                     mwhharmonicflag, mwdharmonicflag, lmcharmonicflag);


   orbit(6,step) = fx;
   orbit(7,step) = fy;
   orbit(8,step) = fz;

   // advance velocities
   for (j=3; j<6; j++) {
     orbit(j,step) = orbit(j,step-1) + (0.5*(orbit(j+3,step-1)+orbit(j+3,step))  * dt );
   }

 }


 // convert to physical units again...
 if (discframe) {
   for (int n=0; n<nint; n++) {
     tvirnow  = tvirbegin - dtvir*n;
     // get the present-day MWD coordinates: the zero of the total
     return_centre(tvirnow, MWD->orient, disccoords);
     //std::cout << setw(14) << tvirnow  << setw(14) << zerocoords[0] << setw(14) << zerocoords[1] << setw(14) << zerocoords[2] << std::endl;
     return_vel_centre(tvirnow, MWD->orient, discvelcoords);
     for (j=0; j<3; j++) {
       orbit(j,n)   = orbit(j,n) - virial_to_physical_length(disccoords[j]);
       orbit(j+3,n) = orbit(j+3,n) - virial_to_physical_velocity(discvelcoords[j]);
     }
   }
 }

 return orbit;

}





std::vector< MatrixXd > MWLMC::rewind(MatrixXd xinit,
                                      MatrixXd vinit,
                                      double dt,
                                      int mwhharmonicflag, int mwdharmonicflag, int lmcharmonicflag,
                                      double rewind_time,
                                      bool discframe,
                                      bool verbose)
{
 /*
 rewind
 -------------
 for a given position and velocity, rewind the orbit by a specified amount of time

 takes physical units

  */

  // step 1: get the coordinate frame of the MW disc
  vector<double> zerocoords(3),zerovelcoords(3);
  return_centre(reference_time, MWD->orient, zerocoords);
  return_vel_centre(reference_time, MWD->orient, zerovelcoords);

  if (verbose) {
    std::cout << setw(14) << reference_time
              << setw(14) << virial_to_physical_length(zerocoords[0])
              << setw(14) << virial_to_physical_length(zerocoords[1])
              << setw(14) << virial_to_physical_length(zerocoords[2]) << std::endl;
    std::cout << setw(14) << "(u,v,w)"
              << setw(14) << virial_to_physical_velocity(zerovelcoords[0])
              << setw(14) << virial_to_physical_velocity(zerovelcoords[1])
              << setw(14) << virial_to_physical_velocity(zerovelcoords[2]) << std::endl;
 }

 double fx,fy,fz,tvir;

 int norbits = xinit.rows();
 std::vector< MatrixXd > orbit;
 orbit.resize(norbits);

 // number of integration steps to take
 int nint = (int)(rewind_time/dt);

 if (nint<0) {
   std::cerr << "rewind: rewind_time must be positive. Exiting." << std::endl;
   exit(-1);
 }


 // set beginning times and timesteps
 double tvirbegin  = reference_time;
 double dtvir      = physical_to_virial_time_return(dt);
 double tphysbegin = 0.;
 double dtphys     = dt;

 // initialise positions and (inverted) velocities for backwards integration
 // include the forces for now
 for (int n=0;n<norbits;n++) {
   orbit[n].resize(10,nint);
   orbit[n](0,0) =  xinit(n,0);
   orbit[n](1,0) =  xinit(n,1);
   orbit[n](2,0) =  xinit(n,2);
   orbit[n](3,0) = -vinit(n,0);
   orbit[n](4,0) = -vinit(n,1);
   orbit[n](5,0) = -vinit(n,2);
   orbit[n](9,0) = tphysbegin;
 }

 //now step forward one, using leapfrog (drift-kick-drift) integrator
 //    https://en.wikipedia.org/wiki/Leapfrog_integration
 //
 int step = 1;

 // get the initial coefficient values: the time here is in tvir units, so always start with 0
 MatrixXd tcoefsmw,tcoefslmc, mwcoscoefs,mwsincoefs;
 MW->select_coefficient_time(tvirbegin, tcoefsmw);
 LMC->select_coefficient_time(tvirbegin, tcoefslmc);
 MWD->select_coefficient_time(tvirbegin, mwcoscoefs, mwsincoefs);

 // return forces for the initial step
 for (int n=0;n<norbits;n++) {
   all_forces_coefs(tcoefsmw, tcoefslmc, mwcoscoefs, mwsincoefs,
                    tphysbegin, orbit[n](0,0),orbit[n](1,0),orbit[n](2,0),
                    fx, fy, fz,
                    mwhharmonicflag, mwdharmonicflag, lmcharmonicflag);

   orbit[n](6,0) = fx;
   orbit[n](7,0) = fy;
   orbit[n](8,0) = fz;
  }

 int j;
 double tvirnow, tphysnow;

 for (step=1; step<nint; step++) {

   // 'advance' timestep: this is in virial units by definition.
   tvirnow  = tvirbegin  - dtvir*step;
   tphysnow = tphysbegin - dtphys*step;

   MW->select_coefficient_time(tvirnow, tcoefsmw);
   LMC->select_coefficient_time(tvirnow, tcoefslmc);
   MWD->select_coefficient_time(tvirnow, mwcoscoefs, mwsincoefs);

   for (int n=0;n<norbits;n++) {
     orbit[n](9,step) = tphysnow; // record the time

     // 'advance' positions
     for (j=0; j<3; j++) {
       orbit[n](j,step) = orbit[n](j,step-1)   + (orbit[n](j+3,step-1)*dt  )  + (0.5*orbit[n](j+6,step-1)  * (dt*dt));
     }

     // this call goes out in physical units
     all_forces_coefs(tcoefsmw, tcoefslmc, mwcoscoefs, mwsincoefs,
                      // need to shift the time by one unit to get the proper centre
                      tphysnow - dtphys, orbit[n](0,step),orbit[n](1,step),orbit[n](2,step),
                      fx, fy, fz,
                      mwhharmonicflag, mwdharmonicflag, lmcharmonicflag);


     orbit[n](6,step) = fx;
     orbit[n](7,step) = fy;
     orbit[n](8,step) = fz;

     // advance velocities
     for (j=3; j<6; j++) {
       orbit[n](j,step) = orbit[n](j,step-1) + (0.5*(orbit[n](j+3,step-1)+orbit[n](j+3,step))  * dt );
     }
   } // orbit loop

 }


 // convert to physical units again...
 // get the original coordinate centre
 vector<double> initcoords(3),initvel(3);
 tvirnow  = tvirbegin;
 return_centre(tvirnow, MWD->orient, initcoords);
 return_vel_centre(tvirnow, MWD->orient, initvel);

 if (discframe) {
   for (int ni=0; ni<nint; ni++) {
     tvirnow  = tvirbegin - dtvir*ni;
     // get the present-day MWD coordinates: the zero of the total
     return_centre(tvirnow, MWD->orient, zerocoords);
     return_vel_centre(tvirnow, MWD->orient, zerovelcoords);
     for (int n=0;n<norbits;n++) {
       for (j=0; j<3; j++) {
         orbit[n](j,ni)   = orbit[n](j,ni) - virial_to_physical_length(zerocoords[j]) + virial_to_physical_length(initcoords[j]);
         orbit[n](j+3,ni) = orbit[n](j+3,ni) - virial_to_physical_velocity(zerovelcoords[j]) + virial_to_physical_length(initvel[j]);
       }
     }
   }
 }

 return orbit;

}


#endif
