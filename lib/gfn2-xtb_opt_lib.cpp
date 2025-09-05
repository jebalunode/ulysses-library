#include "../src/GFN.hpp"
#include <cstdlib>
#include <string>


//program to run GFN2-xTB calculations

std::vector<double> gfn2_xtb_opt(std::string & molecule, int charge=0, double Telec=300.0 ) {
  
  //argument 1 to program is the geometry in xyz format
  //argument 2 to program is the charge
  //argument 3 to program is the electronic temperature (default = 300 K)

  //std::cout << "charge = " << charge << std::endl;
  //std::cout << "T electron = " << Telec << std::endl;
  
  //allocation of molecule with certain charge and multiplicity 2S+1 = 1
  Molecule Mol1(molecule,charge,1);
  
  //allocation of basis set object
  BSet basis(Mol1,"gfn2");
  
  //allocation of GFN2 object
  GFN2 electron(basis,Mol1);
  electron.setElectronTemp(Telec);
  
  //to use ALPB solvation model, uncomment one of the statements below
  electron.setSolvent("water");
  //electron.setSolvent("acetone");
  //electron.setSolvent("acetonitrile");
  //electron.setSolvent("aniline");
  //electron.setSolvent("benzaldehyde");
  //electron.setSolvent("benzene");
  //electron.setSolvent("dichloromethane");
  //electron.setSolvent("chloroform");
  //electron.setSolvent("carbon disulfide");
  //electron.setSolvent("dioxane");
  //electron.setSolvent("dmf");
  //electron.setSolvent("dmso");
  //electron.setSolvent("ethanol");
  //electron.setSolvent("diethyl ether");
  //electron.setSolvent("ethyl acetate");
  //electron.setSolvent("furane");
  //electron.setSolvent("hexadecane");
  //electron.setSolvent("hexane");
  //electron.setSolvent("methanol");
  //electron.setSolvent("nitromethane");
  //electron.setSolvent("octanol");
  //electron.setSolvent("phenol");
  //electron.setSolvent("thf");
  //electron.setSolvent("toluene");
  //electron.setSolvent("water");
  //electron.setSolvent("octanol wet");
  
  //SCF
  electron.Calculate(0);
  BFGSd solve(4,6);

  double energy_threshold = 5.0e-6;
  double gradient_threshold = 1.0e-3;
  SolverOpt(electron,solve,4,0,energy_threshold,gradient_threshold);

  std::cout << std::setprecision(10);

  //get the molecule with optimized geometry
  Molecule Mol2 = electron.Component();


  //get the partial charges
  std::vector<double> AtmCharge = electron.getQAtoms();
  return AtmCharge ;
}
