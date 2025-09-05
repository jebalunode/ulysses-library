#include <iostream>
#include "../src/MNDOd.hpp"
#include "../src/math/SolverPackage.hpp"
#include "../src/Gas.hpp"

int main(int argc, char** argv) {

  //check the files gfn2-xtb_optg.cpp and pm6-corrected.cpp to extend the options here

  //basic stuff
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
  
  //define the method
  BSet basis(Mol1,"pm6");
  PM6 electron(basis,Mol1,"0","D3H4X");
  electron.setEpsilonS(lshift);
  electron.Calculate(0);
  
  //optimize the geometry
  BFGSd solve(4,6);
  SolverOpt(electron,solve,4,0,5e-6,1e-3);
  
  //get vibrations
  std::vector<double> all_vibrations = electron.CalcVibrFrequencies();
  int nvibrations = 6;
  std::vector<double> vibrations(all_vibrations.size() - nvibrations);            //CalcVibrFrequencies returns also translation and rotation modes; these must be removed
  for (size_t idvibr = 0; idvibr < vibrations.size(); ++idvibr) {
    vibrations[idvibr] = all_vibrations[idvibr + nvibrations];
  }
  
  //get electronic energies
  std::vector<double> Eel;
  Eel.push_back(electron.getEnergy(1));      //the one means that the D3H4X correction is applied to the total energy; use 0 if you want non-corrected energies
  //get the degeneracy of ground state
  std::vector<double> gel(1,1.0);
  
  //get the eigenvalues of inertia
  std::vector<double> inertia = electron.Component().InertiaEigenvalues();
  
  double T = 298.15;
  bool grimmecorrection = true;
  double numbermolecules = NA; //1 mol
  double volume = 0.0224;
  PBlRRlHOE IdealGas(T,argv[1],inertia,vibrations,Eel,gel,charge,1,"C1","0",grimmecorrection,numbermolecules,volume);
  
  //loop over temperatures
  double temperature = 100.0;  //K
  std::cout << "temperature (K)      U (J/[K.mol])   " << std::endl;
  for (size_t idx = 0; idx < 100; ++idx) {
    IdealGas.changeT(temperature);
    std::cout << temperature << "  " << IdealGas.S() << std::endl;
    //.U() gets internal energy; replace with
    //.H() -> enthalpy (J/mol)
    //.U() -> internal energy (J/mol)
    //.G() -> Gibbs energy (J/mol)
    //.A() -> Helmholtz energy (J/mol)
    //.P() -> pressure (this is basically PV=nRT; Pa)
    //.CV() -> heat capacity (constant V; J/(K.mol))
    //.CP() -> heat capacity (constant P; J/(K.mol))
    temperature += 10.0;
  }
  
  return 0;
}
