#include "../src/MNDOd.hpp"
#include "../src/math/SolverPackage.hpp"
#include <stdlib.h>

//program to run PM6-D3H4X calculations

int main(int argc, char** argv) {

  //check the files gfn2-xtb_optg.cpp and pm6-corrected.cpp to extend the options here

  std::cout << "running " << argv[1] << "\n";

  char *p;
  int charge = strtol(argv[2],&p,10);
  char *q;
  double lshift = strtod(argv[3],&q);
  
  std::cout << "charge          = " << charge << std::endl;
  std::cout << "level shift     = " << lshift << std::endl;

  Molecule Mol1(argv[1],charge,1,"C1");
  
  std::vector<size_t> atoms = Mol1.Atoms();
  matrixE Geometry = Mol1.Geometry();
  
  BSet basis(Mol1,"pm6");
  
  PM6 electron(basis,Mol1,"0","D3H4X");
  electron.setEpsilonS(lshift);
  
  electron.Calculate(0);
  
  BFGSd solve(4,6);
  SolverOpt(electron,solve,4,0,5e-6,1e-3);
  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  std::string newgeom(argv[1]);
  size_t extpos = newgeom.find(".xyz");
  size_t lengeom = newgeom.size();
  newgeom.replace(extpos,lengeom,"");
  newgeom += "_opt";
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(newgeom);
  
  return 0;
}