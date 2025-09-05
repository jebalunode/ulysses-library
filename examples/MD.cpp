#include <iostream>
#include "../src/MolecularDynamics.hpp"
#include "../src/GFN.hpp"

int main(int argc, char** argv) {
  //Test MD
  std::cout << "Testing Molecular Dynamics -------------------------------" << std::endl;

  //argument 1 to program is the geometry in xyz format
  //argument 2 to program is the charge
  //argument 3 to program is the temperature in Kelvin
  //argument 4 to program is the total simulation time (ps)
  //argument 5 to program is the time time (ps)

  std::cout << "running " << argv[1] << "\n";

  char *p;
  int charge = strtol(argv[2],&p,10);
  char *s;
  double Temperature = strtod(argv[3],&s);
  char *q;
  double tmax = strtod(argv[4],&q);
  char *r;
  double tstep = strtod(argv[5],&r);
  
  double toptg = 0.01;
  bool doequ = false;
  bool optgeometry = true;
  
  Molecule mol(argv[1],charge,1);
  BSet basis(mol,"gfn2");
  GFN2 electron(basis,mol);
  
  Dynamics MDobj(mol,tmax,doequ,optgeometry,tstep,toptg);
  MDobj.ApplyConstraints(false,0);
  MDobj.setIntegration("LeapFrog");
  
  //where to save files
  MDobj.setTrajectoryFile("MD/trajectory");
  MDobj.setEquilibrationFile("MD/equilibration");
  MDobj.setGeometryFile("MD/geometry");
  
  MDobj.runMD(electron,Temperature);
  
  return 0;
}
