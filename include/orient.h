/*
orient.h

functions to
1. read in the orientation files
2. interpolate between different orient setups
3. approximate trajectory outside of simulation bounds (dangerous!)

MSP 22 Apr 2020 clean version
MSP 25 Apr 2020 add linear interpolation for speed
MSP  1 Sep 2020 add spline velocity, warnings
MSP  1 Sep 2020 add new pre-simulation treatment with `backwards' flag
MSP 24 Dec 2021 enable inertial centres if no orient file is passed
MSP 24 Dec 2021 change spline calls to preprocessor flag

wishlist:
-handle bad files gracefully
-handle axis tipping orient files
-find_initial_velocity: how to dynamically determine nPoints?
-offload backwards,backwardsaccel as preprocessor directives

 */
#ifndef ORIENT_H
#define ORIENT_H

#include <iostream>
#include <fstream>

//using namespace std;
using std::cout, std::cerr, std::endl, std::setw;


// decide on splines
// NOTE
// do not use splines unless you are sure it is what you want. the
// spline fits will change the centres in such a way that they are NOT
// guaranteed to match the initial simulation.
#if SPLINEORIENT
#include "spline.h"
#endif


// if true, enable extrapolation before the EXP simulation started
// NOTE
// do NOT use this method unless you know what is going to happen.
bool backwards      = false;
bool backwardsaccel = false;


bool havevelocity   = true; // true if the orient file also has velocities.

struct SphOrient
{

  bool inertial=true;   // stick with inertial centre (default true)?
  bool eventime=true;   // check spacing of orient timing (assume equal)

  int NUMT;             // the number of timesteps
  vector<double> time;  // the time vector, len NUMT
  vector<double> xcen;  // the xcen vector, len NUMT
  vector<double> ycen;  // the ycen vector, len NUMT
  vector<double> zcen;  // the zcen vector, len NUMT

  // velocity centres
  vector<double> ucen;  // the ucen (x-velocity) vector, len NUMT
  vector<double> vcen;  // the vcen (y-velocity) vector, len NUMT
  vector<double> wcen;  // the wcen (z-velocity) vector, len NUMT

#if SPLINEORIENT
  tk::spline xspline;   // spline representation of xcenter
  tk::spline yspline;   // spline representation of ycenter
  tk::spline zspline;   // spline representation of zcenter

  tk::spline vxspline;  // spline representation of vxcenter
  tk::spline vyspline;  // spline representation of vycenter
  tk::spline vzspline;  // spline representation of vzcenter
#endif

  // initial velocity fits
  vector<double> zerotimevelocities;
  vector<double> zerotimeintercepts;

};

void find_time_index(double desired_time, SphOrient orient,
		     int& indx, double& dt)
{
  // starting at the first indx, stop when we get to the matching time
  indx = 0;
  while (orient.time[indx]<=desired_time) {
    indx ++;
  }

  // reset by one
  indx --;

  // guard against wanton extrapolation: should this stop the model?
  if (indx>orient.NUMT-2) cerr << "orient::find_time_index: time after to simulation end selected. setting to latest step." << endl;

  if (indx<0) indx = 0;
  if (indx>orient.NUMT-2) indx = orient.NUMT - 2;

  // check the spacing on orient.time (can be nonuniform)
  dt = orient.time[indx+1] - orient.time[indx];
}


void interpolate_centre(double desired_time,
			SphOrient& orient,
			vector<double>& centre)
{

  // find the matching time
  int indx;
  double dt;

  if (orient.eventime) {
    // the orient array MUST be evenly spaced, check exp output
    dt   = orient.time[1] - orient.time[0];
    indx = (int)( (desired_time-orient.time[0])/dt);
  } else {
    // if the orient array is not evenly spaced, scan until the correct time is found
    find_time_index(desired_time, orient, indx, dt);
  }


  if (havevelocity && indx<0) {
    // interpolate the centre backwards in time before the simulation starts
    // selects the earliest time
    float dtime = (desired_time-orient.time[0]);

    if (backwardsaccel) {
      centre[0] = orient.xcen[0] + dtime*(orient.zerotimevelocities[0]*dtime + orient.zerotimeintercepts[0]);
      centre[1] = orient.ycen[0] + dtime*(orient.zerotimevelocities[1]*dtime + orient.zerotimeintercepts[1]);
      centre[2] = orient.zcen[0] + dtime*(orient.zerotimevelocities[2]*dtime + orient.zerotimeintercepts[2]);
    } else {
      centre[0] = orient.xcen[0] + dtime*orient.zerotimevelocities[0];
      centre[1] = orient.ycen[0] + dtime*orient.zerotimevelocities[1];
      centre[2] = orient.zcen[0] + dtime*orient.zerotimevelocities[2];
    }

  } else {
    // standard mode
    if (indx>orient.NUMT-2) indx = orient.NUMT - 2;

    double x1 = (orient.time[indx+1] - desired_time)/dt;
    double x2 = (desired_time - orient.time[indx])/dt;

    centre[0] = (x1*orient.xcen[indx] + x2*orient.xcen[indx+1]);
    centre[1] = (x1*orient.ycen[indx] + x2*orient.ycen[indx+1]);
    centre[2] = (x1*orient.zcen[indx] + x2*orient.zcen[indx+1]);
  }

}

void interpolate_velocity_centre(double desired_time,
			         SphOrient& orient,
			         vector<double>& velcentre)
{
  // find the matching time
  int indx;
  double dt;

  if (orient.eventime) {
    // the orient array MUST be evenly spaced, check exp output
    dt   = orient.time[1] - orient.time[0];
    indx = (int)( (desired_time-orient.time[0])/dt);
  } else {
    find_time_index(desired_time, orient, indx, dt);
  }

  if (indx<0) {
    // if the desired_time is earlier than the simulation beginning
    // return the approximation for the first velocity

    if (backwardsaccel) {
      float dtime = (desired_time-orient.time[0]);
      velcentre[0] = orient.zerotimevelocities[0]*dtime + orient.zerotimeintercepts[0];
      velcentre[1] = orient.zerotimevelocities[1]*dtime + orient.zerotimeintercepts[1];
      velcentre[2] = orient.zerotimevelocities[2]*dtime + orient.zerotimeintercepts[2];
    } else {
      velcentre[0] = orient.zerotimevelocities[0];
      velcentre[1] = orient.zerotimevelocities[1];
      velcentre[2] = orient.zerotimevelocities[2];
    }

  } else {
    // verify that this works? shouldn't ever get to this point.
    if (indx>orient.NUMT-2) indx = orient.NUMT - 2;

    double x1 = (orient.time[indx+1] - desired_time)/dt;
    double x2 = (desired_time - orient.time[indx])/dt;

    velcentre[0] = (x1*orient.ucen[indx] + x2*orient.ucen[indx+1]);
    velcentre[1] = (x1*orient.vcen[indx] + x2*orient.vcen[indx+1]);
    velcentre[2] = (x1*orient.wcen[indx] + x2*orient.wcen[indx+1]);
  }

}

#if SPLINEORIENT
void spline_centre(double desired_time,
		   SphOrient& orient,
		   vector<double>& centre)
{
  centre[0] = orient.xspline(desired_time);
  centre[1] = orient.yspline(desired_time);
  centre[2] = orient.zspline(desired_time);
}

void spline_vel_centre(double desired_time,
		   SphOrient& orient,
		   vector<double>& velcentre)
{
  // compute the centre of mass velocity, using splines
  velcentre[0] = orient.vxspline(desired_time);
  velcentre[1] = orient.vyspline(desired_time);
  velcentre[2] = orient.vzspline(desired_time);
}
#endif

void return_centre(double desired_time,
		   SphOrient& orient,
		   vector<double>& centre)
{

  if (orient.inertial) {

    // set velocity centre to zero if no orient file
    centre[0] = 0.;
    centre[1] = 0.;
    centre[2] = 0.;

  } else {

#if SPLINEORIENT
      spline_centre(desired_time, orient, centre);
#else
      interpolate_centre(desired_time, orient, centre);
#endif

  }

}

void return_vel_centre(double desired_time,
		       SphOrient& orient,
		       vector<double>& velcentre)
{

  if (orient.inertial) {

    // set velocity centre to zero if no orient file
    velcentre[0] = 0.;
    velcentre[1] = 0.;
    velcentre[2] = 0.;

  } else {
#if SPLINEORIENT
      spline_vel_centre(desired_time, orient, velcentre);
#else
      interpolate_velocity_centre(desired_time, orient, velcentre);
#endif
  }

}

void find_initial_velocity(SphOrient& orient,
			   bool accel=false,
			   int nPoints=2000)
{
  // find the slope and intercept at the beginning of the series, from
  // the first nPoints points

  // use the wikipedia linear regression steps:
  // https://en.wikipedia.org/wiki/Simple_linear_regression#Fitting_the_regression_line

  // and the classic implementation
  // https://stackoverflow.com/questions/11449617/how-to-fit-the-2d-scatter-data-with-a-line-with-c

  double sumT=0, sumX=0, sumY=0, sumZ=0, sumTX=0, sumTY=0, sumTZ=0, sumT2=0;

  vector<double> xterm,yterm,zterm;
  if (accel) {
    xterm = orient.ucen;
    yterm = orient.vcen;
    zterm = orient.wcen;
  } else {
    xterm = orient.xcen;
    yterm = orient.ycen;
    zterm = orient.zcen;
  }

  for(int i=0; i<nPoints; i++) {
      sumT  += orient.time[i];
      sumX  += xterm[i];
      sumY  += yterm[i];
      sumZ  += zterm[i];
      sumTX += orient.time[i] * xterm[i];
      sumTY += orient.time[i] * yterm[i];
      sumTZ += orient.time[i] * zterm[i];
      sumT2 += orient.time[i] * orient.time[i];
  }
  double tMean = sumT / nPoints;
  double xMean = sumX / nPoints;
  double yMean = sumY / nPoints;
  double zMean = sumZ / nPoints;
  double denominator = sumT2 - sumT * tMean;

  if ( std::fabs(denominator) < 1e-7 ) {
    // Fail: it seems a vertical line
    cerr << "orient::find_initial_velocity: can't extrapolate, given these times.";
    // what should the failure look like?
    exit(-1);
  }

  // compute slopes and intercepts
  orient.zerotimevelocities[0]   = (sumTX - sumT * xMean) / denominator;
  orient.zerotimevelocities[1]   = (sumTY - sumT * yMean) / denominator;
  orient.zerotimevelocities[2]   = (sumTZ - sumT * zMean) / denominator;

  if (accel) {
    orient.zerotimeintercepts[0] = xMean - orient.zerotimevelocities[0] * tMean;
    orient.zerotimeintercepts[1] = yMean - orient.zerotimevelocities[1] * tMean;
    orient.zerotimeintercepts[2] = zMean - orient.zerotimevelocities[2] * tMean;
  }

}


void read_orient (string orient_file, SphOrient& orient) {

  // this would be great to modify to take raw exp orient files as well

  if (orient_file=="") {
    cout << "orient::read_orient: no orient file detected . . . proceeding with inertial centre." << endl;
    return;
  }

  ifstream infile;
  infile.open(orient_file, ios::in);
  if (!infile) {
        cout << "orient::read_orient: unable to open specified file. exiting.";
        exit(1); // terminate with error
  }

  orient.inertial = false;

  string line;

  int linenum = 0;

  while (getline(infile, line)) {
    istringstream ss(line);

    if (linenum==0) {
    	ss >> orient.NUMT;
    	orient.time.resize(orient.NUMT);
    	orient.xcen.resize(orient.NUMT);
    	orient.ycen.resize(orient.NUMT);
    	orient.zcen.resize(orient.NUMT);
    	if (havevelocity) {
    	  orient.ucen.resize(orient.NUMT);
    	  orient.vcen.resize(orient.NUMT);
    	  orient.wcen.resize(orient.NUMT);
    	}
    } else {
    	if (havevelocity) {
        ss >>
        orient.time[linenum-1] >>
        orient.xcen[linenum-1] >>
        orient.ycen[linenum-1] >>
    	  orient.zcen[linenum-1] >>
        orient.ucen[linenum-1] >>
        orient.vcen[linenum-1] >>
    	  orient.wcen[linenum-1];
    	} else {
        ss >>
        orient.time[linenum-1] >>
        orient.xcen[linenum-1] >>
        orient.ycen[linenum-1] >>
    	  orient.zcen[linenum-1];
	    }

    }
  linenum ++;

  }

  infile.close();

  // check if time spacing is consistent
  double dt = orient.time[1] - orient.time[0];
  for (int i=2;i<orient.NUMT;i++) {
    if (abs(orient.time[i] - orient.time[i-1] - dt) > dt/10) {
      //cout << setw(14) << i << setw(14) << orient.time[i] << setw(14) << abs(orient.time[i] - orient.time[i-1] - dt) << setw(14) << dt << endl;
      orient.eventime = false;
    }
  }

  // for some future verbose flag, perhaps
#if DEBUGCOEFS
  if (orient.eventime) std::cout << "orient.read_orient: found even time spacing" << "\n";
  std::cout << setw(18) << orient.time[0] << setw(18) << orient.NUMT << "\n";
#endif

#if SPLINEORIENT
    // construct spline interpolations if necessary
    orient.xspline.set_points(orient.time,orient.xcen);
    orient.yspline.set_points(orient.time,orient.ycen);
    orient.zspline.set_points(orient.time,orient.zcen);

    if (havevelocity) {
      // only construct the velocity spline interpolations if REALLY necessary
      orient.vxspline.set_points(orient.time,orient.ucen);
      orient.vyspline.set_points(orient.time,orient.vcen);
      orient.vzspline.set_points(orient.time,orient.wcen);
    }
#endif

  if (havevelocity) {
    orient.zerotimevelocities.resize(3);

    if (backwards) {
      // enable extrapolation to before the simulation started
      // (but see warning at top of file)

      if (backwardsaccel) {
        orient.zerotimeintercepts.resize(3);
        find_initial_velocity(orient,true);
      } else {
        find_initial_velocity(orient);
      }

    } else {
      // set the initial velocities to be the t=0 velocities
      orient.zerotimevelocities[0] = orient.ucen[0];
      orient.zerotimevelocities[1] = orient.vcen[0];
      orient.zerotimevelocities[2] = orient.wcen[0];
    }
  } // end havevelocity

}

#endif