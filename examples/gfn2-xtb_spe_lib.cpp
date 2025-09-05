#include "../src/GFN.hpp"
#include <cstdlib>
#include <string>


//program to run GFN2-xTB calculations

std::vector<double>  gfn2_xtb_spe(std::string & molecule, int charge, double Telec ) {
  
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
  //get the partial charges
  std::vector<double> AtmCharge = electron.getQAtoms();
  return AtmCharge ;
  //get the atom list
  //std::vector<size_t> atoms = Mol1.Atoms();
  //size_t Natoms = atoms.size();
  //AO basis info
  //std::vector<size_t> AOS = basis.AtomNAOs(atoms);
  
  //get the atomic polarizabilities
  //std::vector<double> polarizabilities;
  // electron.AtomicPolarizabilities(polarizabilities,AtmCharge);
  
  //print
  /*
  std::cout << std::setprecision(5) << "\n";
  //std::cout << "atom       AOs          charge           pol\n";
  std::cout<<"charge\n";
  for (size_t idx = 0; idx < Natoms; ++idx) {
    //std::cout << AtomNr2Symbol(atoms[idx]) << "          ";
    //std::cout << AOS[idx] << "          ";
  //  if (AtmCharge[idx] > 0.0) {std::cout << " ";}
    std::cout << AtmCharge[idx] <<"\n";
    //<< "          " << polarizabilities[idx] << "\n";
  }
  std::cout << "\n";
  */
}
