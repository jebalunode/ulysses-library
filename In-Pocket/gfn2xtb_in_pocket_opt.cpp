#include "../src/GFN.hpp"
#include "../src/MolecularDynamics.hpp"
#include "../src/math/SolverPackage.hpp"
#include <stdlib.h>

//program to run GFN2-xTB calculations with in-pocket

int main(int argc, char** argv) {
  
  //arguments
  // 0 exe
  // 1 geometry
  // 2 charge
  // 3 alpha = 1.0
  // 4 kappa = 0.025
  // 5 ignore protons?
  // 6 use unconstrained optimisation?
  // 7 Telec
  // 8 solvation?
  // 9 solvent name
  // 10 calculate density?
  // 11 name of density file
  // 12 name of geometry file
  // 13 energy threshold for geometry optimization = 1.0e-6
  // 14 gradient threshold for geometry optimization = 1.0e-3
  
  //parameters passed as argument
  char *p;
  int charge = strtol(argv[2],&p,10);
  char *r;
  double alpha = strtod(argv[3],&r);
  char *s;
  double kappa = strtod(argv[4],&s);
  char *t;
  int ignh = strtol(argv[5],&t,10);
  char *u;
  int constrainedopt = strtol(argv[6],&u,10);
  char *q;
  double Telec = strtod(argv[7],&q);
  char *v;
  int solvation = strtol(argv[8],&v,10);
  char *w;
  int calcdensity = strtol(argv[10],&w,10);
  char *y;
  double energy_threshold = strtod(argv[13],&y);
  char *z;
  double gradient_threshold = strtod(argv[14],&z);
  
  int natomspotential = 0;

  std::cout << "running " << argv[1] << "\n";
  std::cout << "charge     = " << charge << std::endl;
  std::cout << "alpha      = " << alpha << std::endl;
  std::cout << "kappa      = " << kappa << std::endl;

  Molecule Mol1(argv[1],charge,1,"C1");

  std::vector<size_t> atoms = Mol1.Atoms();
  size_t Natoms = atoms.size();
  matrixE Geometry = Mol1.Geometry();

  BSet basis(Mol1,"gfn2");
  std::vector<size_t> AOS = basis.AtomNAOs(atoms);

  GFN2 electron(basis,Mol1);
  electron.setElectronTemp(Telec);
  
  //use ALPB solvation?
  if (solvation > 0) {
    electron.setSolvent(argv[9]);
    // "water"
    // "acetone"
    // "acetonitrile"
    // "aniline"
    // "benzaldehyde"
    // "benzene"
    // "dichloromethane"
    // "chloroform"
    // "carbon disulfide"
    // "dioxane"
    // "dmf"
    // "dmso"
    // "ethanol"
    // "diethyl ether"
    // "ethyl acetate"
    // "furane"
    // "hexadecane"
    // "hexane"
    // "methanol"
    // "nitromethane"
    // "octanol"
    // "phenol"
    // "thf"
    // "toluene"
    // "water"
    // "octanol wet"
  }

  electron.setRestart(0);                                      //use restart, otherwise this will be a killer
  electron.Calculate(0);

  BFGSd solve(4,6);
  MetaDynamics MTD;
  if (constrainedopt == 0) {
    SolverOpt(electron,solve,4,0,energy_threshold,gradient_threshold);
  }
  else {
    MTD.setMetaSet2Structure(Geometry);
    MTD.setMetaAlphas2Alpha(alpha);
    for (size_t idAtm = 0; idAtm < Natoms; ++idAtm) {
      if ((ignh == 0)&&(atoms[idAtm] == 1)) {continue;}
      MTD.addMetaAtom(idAtm + 1);
      ++natomspotential;
    }
    MTD.setMetaFactors2kPull(-fabs(kappa)*double(natomspotential));       //make sure it is negative
    std::cout << "fixed atoms" << std::endl;
    std::vector<size_t> metatomsprint = MTD.MetaAtoms();
    for (size_t idAtm = 0; idAtm < metatomsprint.size(); ++idAtm) {
      std::cout << metatomsprint[idAtm] << " ";
    }
    std::cout << std::endl;
  
    CompositeMTDQM mixmethod(MTD,electron,Geometry);             //this contains pointers to the two methods, so the object should be always up-to-date
    SolverOpt(mixmethod,solve,4,0,energy_threshold,gradient_threshold);
  }
  std::cout << std::setprecision(10);
  
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(argv[12]);
  
  electron.Calculate(0);
  
  matrixE optGeometry = electron.Component().Geometry();
  
  std::cout << "Optimized g energy: " << electron.getEnergy() << " Hartree\n\n";
  
  matrixE aux = optGeometry;
  double rmsd = MTD.RMSD(Geometry,optGeometry,aux);
  
  std::cout << "RMSD with initial structure: " << rmsd << " \n\n";
  
  std::vector<double> AtmCharge = electron.getQAtoms();
  
  std::vector<double> polarizabilities;
  electron.AtomicPolarizabilities(polarizabilities,AtmCharge);
  
  std::cout << std::setprecision(5) << "\n";
  std::cout << ">atom;charge;pol\n";
  for (size_t idx = 0; idx < Natoms; ++idx) {
    std::cout << atoms[idx] << "          ";
    if (atoms[idx] < 10) {std::cout << " ";}
    std::cout << AOS[idx] << "          ";
    if (AtmCharge[idx] > 0.0) {std::cout << " ";}
    std::cout << AtmCharge[idx] << "          " << polarizabilities[idx] << "\n";
  }
  std::cout << "<atom;charge;pol\n";
  
  double polbity = 0.0;
  electron.TotalPolarizability(polbity,AtmCharge);
  std::cout << " Total Polarizability          " << polbity << "\n";
  
  //get the density?
  if (calcdensity > 0) {
    electron.ElectronicDensity(argv[11]);
  }
  
  return 0;
}