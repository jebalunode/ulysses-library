#include <iostream>
#include "../src/GFN.hpp"
#include "../src/math/SolverPackage.hpp"
#include "../src/Gas.hpp"

int main(int argc, char** argv) {

  //check the files gfn2-xtb_optg.cpp and pm6-corrected.cpp to extend the options here

  //arguments
  // 0 exe
  // 1 geometry file
  // 2 charge
  // 3 Telec = 300
  // 4 solvation?
  // 5 solvent name
  // 6 optimise geometry?
  // 7 name for new geometry file
  // 8 thermo?
  // 9 energy threshold for geometry optimization = 1.0e-6
  // 10 gradient threshold for geometry optimization = 1.0e-3
  // 11 calculate density?
  // 12 name for density file
  // 13 electronic reactivity indices?
  // 14 orbital reactivity indices?
  // 15 Koopman IP?
  // 16 IP?
  // 17 EA?
  // 18 electronativity?
  // 19 hardness?
  
  //parameters passed as argument
  char *p;
  int charge = strtol(argv[2],&p,10);
  char *q;
  double Telec = strtod(argv[3],&q);
  char *r;
  int solvation = strtol(argv[4],&r,10);
  char *v;
  int optgeom = strtol(argv[6],&v,10);
  char *s;
  int thermo = strtol(argv[8],&s,10);
  char *t;
  double energy_threshold = strtod(argv[9],&t);
  char *u;
  double gradient_threshold = strtod(argv[10],&u);
  char *w;
  int calcdensity = strtol(argv[11],&w,10);
  char *z1;
  int elecrx = strtol(argv[13],&z1,10);
  char *z2;
  int orbrx = strtol(argv[14],&z2,10);
  char *z3;
  int koopman = strtol(argv[15],&z3,10);
  char *z4;
  int ip = strtol(argv[16],&z4,10);
  char *z5;
  int ea = strtol(argv[17],&z5,10);
  char *z6;
  int electronegativity = strtol(argv[18],&z6,10);
  char *z7;
  int hardness = strtol(argv[19],&z7,10);
  
  //system declaration
  std::cout << "running " << argv[1] << "\n";
  std::cout << "charge          = " << charge << std::endl;
  std::cout << "T electron = " << Telec << std::endl;

  //allocate molecules
  Molecule Mol1(argv[1],charge,1,"C1");
  
  //define method and basis set
  BSet basis(Mol1,"gfn2");
  GFN2 electron(basis,Mol1);
  electron.setElectronTemp(Telec);
  
  //use ALPB solvation?
  if (solvation > 0) {
    electron.setSolvent(argv[5]);
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

  electron.Calculate(0);
  
  //optimise geometry?
  if (optgeom > 0) {
    BFGSd solve(4,6);
    SolverOpt(electron,solve,4,0,energy_threshold,gradient_threshold);
  }
  Molecule Mol2 = electron.Component();
  Mol2.WriteXYZ(argv[7]);
  
  electron.Calculate(1);
  
  std::cout << std::setprecision(7) << "\n";
  //perform thermodynamics?
  if (thermo > 0) {
    //get vibrations
    std::vector<double> all_vibrations = electron.CalcVibrFrequencies();
    int nvibrations = 6;
    std::vector<double> vibrations(all_vibrations.size() - nvibrations);            //CalcVibrFrequencies returns also translation and rotation modes; these must be removed
    for (size_t idvibr = 0; idvibr < vibrations.size(); ++idvibr) {
      vibrations[idvibr] = all_vibrations[idvibr + nvibrations];
    }
    std::cout << ">all vibrational frequencies" << std::endl;
    for (size_t idvibr = 0; idvibr < all_vibrations.size(); ++idvibr) {
      std::cout << all_vibrations[idvibr] << std::endl;
    }
    std::cout << "<all vibrational frequencies" << std::endl;
    //get electronic energies
    std::vector<double> Eel;
    Eel.push_back(electron.getEnergy(1));      //the one means that the D3H4X correction is applied to the total energy; use 0 if you want non-corrected energies
    //get the degeneracy of ground state
    std::vector<double> gel(1,1.0);
    
    //get the eigenvalues of inertia matrix
    std::vector<double> inertia = electron.Component().InertiaEigenvalues();
    
    double T = 298.15;
    bool grimmecorrection = true;
    double numbermolecules = NA; //1 mol
    double volume = 0.0224;
    PBlRRlHOE IdealGas(T,argv[1],inertia,vibrations,Eel,gel,charge,1,"C1","0",grimmecorrection,numbermolecules,volume);
    
    //loop over temperatures and print out
    double temperature = 100.0;  //K
    std::cout << ">Thermodynamics" << std::endl;
    for (size_t idx = 0; idx < 2201; ++idx) {
      IdealGas.changeT(temperature);
      std::cout << temperature << ";" << IdealGas.S() << ";" << IdealGas.H() << ";" << IdealGas.G() << ";" << IdealGas.U() << ";" << IdealGas.A() << ";" << IdealGas.CP() << ";" << IdealGas.CV() << std::endl;
      temperature += 0.5;
    }
    std::cout << "<Thermodynamics" << std::endl;
  }
  
  //get the density?
  if (calcdensity > 0) {
    electron.ElectronicDensity(argv[12]);
  }

  //get charges and polarisabilities
  std::vector<size_t> atoms = Mol1.Atoms();
  std::vector<double> AtmCharge = electron.getQAtoms();
  std::vector<double> polarizabilities;
  electron.AtomicPolarizabilities(polarizabilities,AtmCharge);
  std::cout << ">atom;charge;pol\n";
  for (size_t idx = 0; idx < atoms.size(); ++idx) {
    std::cout << atoms[idx] << ";";
    std::cout << AtmCharge[idx] << ";" << polarizabilities[idx] << "\n";
  }
  std::cout << "<atom;charge;pol\n";
  double polbity = 0.0;
  electron.TotalPolarizability(polbity,AtmCharge);
  std::cout << " Total Polarizability          " << polbity << "\n";

  //additional properties
  matrixE RxData(1,1);
  if (elecrx > 0) {
    electron.ReactivityIndices(RxData,false);
    std::cout << ">Electronic Reactivity indices" << std::endl;
    RxData.Print(4);
    std::cout << "<Electronic Reactivity indices" << std::endl;
  }

  if (orbrx > 0) {
    electron.ReactivityIndices(RxData,true);
    std::cout << ">Orbital Reactivity indices" << std::endl;
    RxData.Print(4);
    std::cout << "<Orbital Reactivity indices" << std::endl;
  }

  if (koopman > 0) {std::cout << "Ionization Potential (Koopman): " << electron.IonizationPotential(true)*au2eV << "   eV" << std::endl;}
  if (ip > 0) {std::cout << "Ionization Potential (Definition): " << electron.IonizationPotential(false)*au2eV << "   eV" << std::endl;}
  if (ea > 0) {std::cout << "Electron Affinity (Definition): " << electron.ElectronAffinity()*au2eV << "   eV" << std::endl;}
  
  if ((electronegativity > 0)||(hardness > 0)) {
    double chi;
    double eta;
    electron.HSABdata(chi,eta);
    std::cout << "Electronegativity: " << chi*au2eV << "   eV" << std::endl;
    std::cout << "Hardness: " << eta*au2eV << "   eV" << std::endl;
  }
  
  return 0;
}
