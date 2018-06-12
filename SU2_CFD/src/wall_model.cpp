/*!
 * \file wall_model.cpp
 * \brief Function for the wall model functions for hom large eddy simulations.
 * \author E. van der Weide, T. Economon, P. Urbanczyk
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/wall_model.hpp"
#include <lapacke.h>
//#include "../include/solver_structure.hpp"

CWallModel::CWallModel(void){
  thickness = 0.0;
}

CWallModel::~CWallModel(void){}

void CWallModel::Initialize(CBoundaryFEM * boundary, CConfig *config, CGeometry *geometry){}

void CWallModel::SetUpExchange(CBoundaryFEM * boundary, CConfig *config, CGeometry *geometry){}

void CWallModel::SetPoints(CSurfaceElementFEM * thisElem, vector<su2double> exchangeCoords){}

void CWallModel::SolveCoupledSystem(std::vector<su2double> exVals, CSurfaceElementFEM * curFace, unsigned short nDim){}

CWallModel1DEQ::CWallModel1DEQ(void) : CWallModel(){
  expansionRatio = 0.0;
  numPoints = 0.0;
  isIsoThermalWall = false;
  isHeatFluxWall = false;
  specWallTemp = 0.0;
  specWallHeatFlux = 0.0;
}

CWallModel1DEQ::~CWallModel1DEQ(void){}

void CWallModel1DEQ::Initialize(CBoundaryFEM * boundary, CConfig * config, CGeometry * geometry){

  /*--- Set up the wall model for this boundary marker ---*/

  /*--- Get the marker name for this boundary ---*/
  std::string markerName = boundary->markerTag;

  /*--- To Do: Determine what kind of wall this is (isothermal or heat flux)
   * If this is an isothermal wall, set isIsoThermalWall to true and
   * get the specified wall temperature. Otherwise, set isHeatFluxWall to true
   * and get the specified wall heatflux. ---*/

  // FOR NOW, ASSUME ISOTHERMAL WALL
  isIsoThermalWall = true;
  specWallTemp = config->GetIsothermal_Temperature(markerName);

  /*--- Get double data from the configuration class object. This should include
   * the wall model thickness and the expansion ratio of the points ---*/
  su2double * tempDoubleData = config->GetWallFunction_DoubleInfo(markerName);
  thickness = tempDoubleData[0];
  expansionRatio = tempDoubleData[1];

  initCondExists = false;

  // Output values for debugging purposes
  std::cout << "Thickness = " << thickness << std::endl;
  std::cout << "Expansion Ratio = " << expansionRatio << std::endl;

  /*--- Get the integer data from the configuration class object. This includes
   * the number of points to be used in the wall model ---*/
  unsigned short int * tempIntData = config->GetWallFunction_IntInfo(markerName);
  numPoints = tempIntData[0];

  // Output the number of points for debugging purposes
  std::cout << "Number of WM Points = " << numPoints << std::endl;
  std::cout << "*******************" << std::endl;

  /*--- Initialize some things on the boundary, itself ---*/
  boundary->wallModelBoundary = true;
  boundary->nWallModelPoints = numPoints;
  boundary->wallModelExpansionRatio = expansionRatio;
  boundary->wallModelThickness = thickness;

  this->SetUpExchange(boundary, config, geometry);

}

void CWallModel1DEQ::SetUpExchange(CBoundaryFEM * boundary, CConfig * config, CGeometry * geometry){

  /*--- This method will set up the exchange locations in each boundary cell. ---*/

  /*--- Create the DGGeometry pointer by casting the geometry as a CMeshFEM_DG ---*/
  CMeshFEM_DG *DGGeometry = dynamic_cast<CMeshFEM_DG *>(geometry);

  // Get number of dimensions
  unsigned short nDim = DGGeometry->GetnDim();

  // Get array of volume elements
  CVolumeElementFEM * volElements = DGGeometry->GetVolElem();

  /*--- Loop over the boundary elements ---*/
  unsigned short numElems = boundary->surfElem.size();
  for(unsigned short iSurfElem = 0; iSurfElem < numElems; ++iSurfElem) {

    // Get a pointer to current surface element
    CSurfaceElementFEM * curElem = &(boundary->surfElem[iSurfElem]);

    // Get the number of degrees of freedom of this face and the corresponding element
    const unsigned short nDofsFace = curElem->DOFsSolFace.size();
    const unsigned short nDofsElement = curElem->DOFsSolElement.size();
    const unsigned long volID = curElem->volElemID;

    // Set the size of the exchange point ID vector and initialize to zero
    curElem->exchangePointIDs.assign(nDofsFace,0);

    // Get the current face's corresponding volume element
    CVolumeElementFEM * curVolume = &(volElements[volID]);

    // Get the current corresponding volume element's polynomial degree
    unsigned short nPoly = curVolume->nPolySol;

    // Step through the dofs of the surface element
    for(unsigned short iDof = 0; iDof < nDofsFace; ++iDof){

      // Set the initial difference to be the specified thickness
      su2double difference = this->thickness;

      // For each surface DOF, there will be nPoly solution DOFs "above" it in the wall-normal direction
      for(unsigned short i = 0; i < nPoly; ++i){
        // We want to look at the nPoly solution DOFs above the current surface DOF. These solution DOFs will be
        // at the current index plus (nPoly+1)^2, 2*(nPoly+1)^2, 3*(nPoly+1)^2, etc in the vector.

        // Calculate the index of the ith solution DOF above the current surface DOF
        unsigned short thisSolDOFIndex = iDof + (i+1)*((nPoly+1)*(nPoly+1));

        // Get the wall distance of this solution DOF
        su2double thisWallDist = curVolume->wallDistanceSolDOFs[thisSolDOFIndex];

        // Compare the wall distance of this DOF with that specified by the wall model for this wall
        // If the difference is smaller than the current difference, this point is closer to the specified
        // wall thickness
        if( fabs(this->thickness - thisWallDist) < difference ){
          // Set the difference between the specified wall model thickness and this solution DOF to be the new difference
          difference = fabs(this->thickness - thisWallDist);

          // Set the current surface DOF's exchange point ID to that of the current solution DOF
          curElem->exchangePointIDs[iDof] = thisSolDOFIndex;
        }
      }
    }
    su2double exchangeHeight = curVolume->wallDistanceSolDOFs[curElem->exchangePointIDs[0]];
    this->SetPoints(curElem,exchangeHeight);
  }
}

void CWallModel1DEQ::SetPoints(CSurfaceElementFEM * thisElem, su2double exchangeHeight){
  /*--- Allocate memory/vector of wall model points for this element. This assumes that there is only
   * one vector of wall model points needed per cell. This assumption is based on the volume cell being
   * a regular hex or triangular prism, where the distances from each wall integration point to its
   * corresponding exchange locations is the same throughout this boundary cell. This may not be a good
   * assumption in general cases, but will work for plane channel flow for now. ---*/

  thisElem->wallModelPoints.resize(this->numPoints);

  /*--- for now, we are just using the y-coordinate ---*/
  /*--- Use a geometric expansion to get wall model point coords between the wall and the exchange location ---*/
  /*--- Find the first cell thickness ---*/
  su2double firstCellThickness = 0.0;
  if(this->expansionRatio == 1.0){
    firstCellThickness = exchangeHeight/(this->numPoints-1);
  }
  else{
    firstCellThickness = exchangeHeight * (1 - this->expansionRatio) / (1 - std::pow(this->expansionRatio,this->numPoints-1));
  }
//  // Output for debugging
//  if(thisElem->boundElemIDGlobal == 0)
//  {
//    std::cout << "For boundElemIDGlobal = " << thisElem->boundElemIDGlobal << std::endl;
//    std::cout << "Specified WM thickness = " << this->thickness << std::endl;
//    std::cout << "Actual exchange wall distance = " << exchangeHeight << std::endl;
//  }

  // Set first wall model point to be at the wall (y=0)
  thisElem->wallModelPoints[0] = 0.0;
  if(this->expansionRatio == 1.0){
    for(unsigned short iPoint = 1; iPoint < numPoints; ++iPoint){
      thisElem->wallModelPoints[iPoint] = thisElem->wallModelPoints[iPoint-1] + exchangeHeight/(this->numPoints-1);
    }
  }
  else{
    for(unsigned short iPoint = 1; iPoint < numPoints; ++iPoint){
      thisElem->wallModelPoints[iPoint] = thisElem->wallModelPoints[iPoint-1] + firstCellThickness * std::pow(this->expansionRatio,iPoint-1);
    }
  }
  // Output for debugging
  if(thisElem->boundElemIDGlobal == 0){
    for(unsigned short i = 0; i < numPoints; ++i){
      std::cout << "y[" << i << "] = " << thisElem->wallModelPoints[i] << std::endl;
    }
  }
}

void CWallModel1DEQ::SolveCoupledSystem(std::vector<su2double> exVals, CSurfaceElementFEM * curFace, unsigned short nDim){
  /*--- This function solves the coupled systems ---*/

  /*--- First, set up vectors and initialize---*/
  //*********************TO DO***************//
  // Move these to be class members
  // Reset them for each new point
  std::vector<su2double> y(numPoints,0.0);
  std::vector<su2double> u(numPoints,0.0);
  std::vector<su2double> T(numPoints,0.0);
  std::vector<su2double> rho(numPoints,0.0);
  std::vector<su2double> mu(numPoints,0.0);
  std::vector<su2double> muTurb(numPoints,0.0);

  su2double tauWall = 0.0;
  su2double rho_ex = 0.0;
  su2double rho_wall = 0.0;
  su2double u_ex = 0.0;
  su2double v_ex = 0.0;
  su2double w_ex = 0.0;
  su2double e_ex = 0.0;

  // Set some constants, assuming air at standard conditions
  su2double C_1 = 1.458e-6;
  su2double S = 110.4;
  su2double R = 287.058;
  su2double kappa = 0.41;
  su2double A = 17;
  su2double gamma = 1.4;
  su2double Pr_lam = 0.7;
  su2double Pr_turb = 0.9;

  bool converged = false;
  bool initCondExists = false;

  // Calculate additional exchange location values
  su2double c_p = (gamma*R)/(gamma-1);
  su2double c_v = R/(gamma-1);
  if(nDim == 2)
  {
    rho_ex = exVals[0];
    u_ex = exVals[1];
    v_ex = exVals[2];
    e_ex = exVals[3];
  }
  else if(nDim == 3)
  {
    rho_ex = exVals[0];
    u_ex = exVals[1];
    v_ex = exVals[2];
    w_ex = exVals[3];
    e_ex = exVals[4];
  }

  // Get velocity magnitude and direction
  su2double veloc_mag_ex = std::sqrt(u_ex*u_ex + w_ex*w_ex);
  su2double veloc_dir_ex = std::acos(u_ex / veloc_mag_ex);

  // Assume calorically perfect gas.
  su2double T_ex = e_ex/c_v;
  su2double P_ex = rho_ex * R *T_ex;

  // Output for debugging
  if(curFace->volElemID == 0){
    std::cout << "i = ";
    for(unsigned short i = 0; i<numPoints; i++){
      std::cout << i << ", ";
    }
    std::cout << std::endl;
  }

  unsigned short j = 0;
  while( converged == false ){

    /*--- Set initial condition if it doesn't already exist ---*/
    if( initCondExists == false ){
      if(curFace->volElemID == 0){
        std::cout << "Setting Initial Condition for WM" << std::endl;
        std::cout << "veloc_mag = " << veloc_mag_ex << std::endl;
        std::cout << "T_ex = " << T_ex << std::endl;
      }
      // Set the initial friction length
      su2double mu_init = C_1 * std::pow(T_ex,1.5) / (T_ex + S);
      tauWall = 0.1 * mu_init * veloc_mag_ex / curFace->wallModelPoints[1];
      su2double u_tau = std::sqrt(tauWall/rho_ex);
      su2double nu_init = mu_init / rho_ex;
      su2double l_tau = nu_init / u_tau;
      for(unsigned short i = 0; i<numPoints; i++){
        // Initial condition is uniform temperature
        T[i] = T_ex;

        /*--- Set the viscosity (constant to begin) ---*/
        mu[i] = mu_init;

        // Get the y-coordinate of the 1-d points from the current face
        y[i] = curFace->wallModelPoints[i];

        /*--- Set the turbulent viscosity ---*/
        su2double D = std::pow(1-std::exp((-y[i]/l_tau)/A),2.0);
        muTurb[i] = kappa * rho_ex * y[i] * u_tau * D;
        //muTurb[i] = 0.0;
      }
      initCondExists = true;
    }

    // Set up momentum equation matrix and rhs(interior points)
    std::vector<su2double> lower(numPoints-1,0.0);
    std::vector<su2double> upper(numPoints-1,0.0);
    std::vector<su2double> diagonal(numPoints,0.0);
    std::vector<su2double> rhs(numPoints,0.0);

    /*--- Solve the differential equation
     * d/dy[ (mu + mu_turb) * du/dy) = 0
     * Re-write as
     * d^2u/dy^2 + g(y)*du/dy = 0
     * where g(y) = f'(y) / f(y)
     * f(y) = (mu + mu_turb) ---*/

    for(unsigned short i=1; i<numPoints-1; i++)
    {
      // For equally spaced points, these are the same
      su2double dy_minus = (y[i] - y[i-1]);
      su2double dy_plus = (y[i+1] - y[i]);

      su2double f_minus = (mu[i-1] + muTurb[i-1]);
      su2double f = (mu[i] + muTurb[i]);
      su2double f_plus = (mu[i+1] + muTurb[i+1]);

      su2double f_prime = (f_plus - f_minus) / (dy_minus + dy_plus);

      su2double g = f_prime/f;

      //lower[i] = (1/(dy_minus*dy_plus) - g/(dy_minus+dy_plus));
      //upper[i] = (1/(dy_plus*dy_plus) + g/(dy_minus+dy_plus));
      //diagonal[i] = -2.0/(dy_minus*dy_plus);
      lower[i] = (1 - 0.5*dy_minus*g);
      upper[i] = (1 + 0.5*dy_minus*g);
      diagonal[i] = -2.0;
      rhs[i] = 0.0;
    }

    /*--- Set boundary conditions in matrix and rhs ---*/
    diagonal[0] = 1.0;
    rhs[0] = 0.0;
    diagonal[numPoints-1] = 1.0;
    rhs[numPoints-1] = veloc_mag_ex;

    // Solve the matrix problem to get the velocity field
    //********LAPACK CALL*******
    int info = 0;
    info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,numPoints,1,lower.data(),diagonal.data(),upper.data(),rhs.data(),numPoints);

    u = rhs;

    if(curFace->volElemID == 0)
    {
      std::cout << "info = " << info << std::endl;
      std::cout << "u = ";
      for(unsigned short i=0; i<numPoints; i++){
        std::cout << u[i] << ",  ";
      }
      std::cout << std::endl;
      std::cout << "mu = ";
      for(unsigned short i=0; i<numPoints; i++){
        std::cout << mu[i] << ",  ";
      }
      std::cout << std::endl;
      std::cout << "muTurb = ";
      for(unsigned short i=0; i<numPoints; i++){
        std::cout << muTurb[i] << ",  ";
      }
      std::cout << std::endl;
    }

    // Set up energy equation matrix and rhs
    for(unsigned short i=1; i<numPoints-1; i++)
    {
      su2double dy_minus = (y[i] - y[i-1]);
      su2double dy_plus = (y[i+1] - y[i]);

      su2double f_minus = c_p * (mu[i-1]/Pr_lam + muTurb[i-1]/Pr_turb);
      su2double f = c_p * (mu[i]/Pr_lam + muTurb[i]/Pr_turb);
      su2double f_plus = c_p * (mu[i+1]/Pr_lam + muTurb[i+1]/Pr_turb);

      su2double f_prime = (f_plus - f_minus) / (dy_minus + dy_plus);

      su2double g = f_prime/f;

      //lower[i] = (1/dy_minus - 0.5*g);
      //upper[i] = (1/dy_plus + 0.5*g);
      //diagonal[i] = -2.0/dy_minus;

      lower[i] = (1 - 0.5*dy_minus*g);
      upper[i] = (1 + 0.5*dy_minus*g);
      diagonal[i] = -2.0;

      su2double rhs_minus = (mu[i-1] + muTurb[i-1]) * u[i-1] * (u[i] - u[i-1]) / dy_minus;
      su2double rhs_plus = (mu[i+1] + muTurb[i+1]) * u[i+1] * (u[i+1] - u[i]) / dy_plus;
      rhs[i] = -1.0*(rhs_plus - rhs_minus)/(dy_minus + dy_plus);
    }

    // Set up energy boundary conditions
    if(isIsoThermalWall == true)
    {
      diagonal[0] = 1.0;
      rhs[0] = specWallTemp;
      diagonal[numPoints-1] = 1.0;
      rhs[numPoints-1] = T_ex;
    }
    else if(isHeatFluxWall == true)
    {
      // TO DO: Set BCs for heat flux condition
    }

    // Solve the matrix problem to get the temperature field
    // *******LAPACK CALL********
    info = 0;
    info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,numPoints,1,lower.data(),diagonal.data(),upper.data(),rhs.data(),numPoints);

    T = rhs;

    if(curFace->volElemID == 0)
    {
      std::cout << "info = " << info << std::endl;
      std::cout << "T = ";
      for(unsigned short i=0; i<numPoints; i++){
        std::cout << T[i] << ",  ";
      }
      std::cout << std::endl;
    }

    su2double tauWall = mu[1] * u[1] / y[1];
    su2double heatFlux = 0.0;
    if(isIsoThermalWall == true){
      heatFlux = c_p*(mu[1]/Pr_lam) * (T[1] - specWallTemp) / y[1];
    }

    if(curFace->volElemID == 0)
    {
      std::cout << "tauWall = " << tauWall << std::endl;
      std::cout << "heatFlux = " << heatFlux << std::endl;
      std::cout << "---------" << std::endl;
    }

    // Given new wall shear and temperature profile, update solution
    rho_wall = P_ex / (R * T[0]);
    su2double u_tau = std::sqrt(tauWall/rho_wall);
    su2double mu_wall = C_1 * std::pow(T[0],1.5) / (T[0] + S);
    su2double nu_wall = mu_wall / rho_wall;
    su2double l_tau = nu_wall / u_tau;
    for(unsigned short i=0; i<numPoints; i++){
      // Set molecular viscosity based on new temp profile
      mu[i] = C_1 * std::pow(T[i],1.5) / (T[i] + S);

      // Set the density based on the new temp profile
      rho[i] = P_ex / (R*T[i]);
      /*--- Set the turbulent viscosity ---*/
      su2double D = std::pow(1-std::exp((-y[i]/l_tau)/A),2.0);
      muTurb[i] = kappa * rho[i] * y[i] * u_tau * D;
      //muTurb[i] = 0.0;
    }

    if(j == 5){
      converged = true;
    }
    j++;
  }
}

