#include "../src/MNDOd.hpp"
#include <stdlib.h>

//program to run PM6-D3H4X calculations

int main(int argc, char** argv) {
  
  //first argument is the xyz file of molecule
  std::cout << "running " << argv[1] << "\n";

  char *p;
  //second argument is charge
  int charge = strtol(argv[2],&p,10);
  char *q;
  //third argument is level shift, use 0.1 as default
  double lshift = strtod(argv[3],&q);
  
  std::cout << "charge          = " << charge << std::endl;
  std::cout << "level shift     = " << lshift << std::endl;
  
  //allocate molecule, singlet
  Molecule Mol1(argv[1],charge,1);
  
  std::vector<size_t> atoms = Mol1.Atoms();
  matrixE Geometry = Mol1.Geometry();
  
  //allocate basis set object
  //replace "pm6" with:
  //"am1", "pm3", "mndo", "mndod",...
  BSet basis(Mol1,"pm6");
  
  //allocate PM6 object with the D3H+ correction
  //replace "D3H+" with "D3H4X" for D3H4X; omit for no correction
  //replace "0" with "UHF" for radicals
  //To get other methods, replace PM6 with AM1, PM3, MNDO, MNDOd, ...
  PM6 electron(basis,Mol1,"0","D3H+");
  
  electron.setEpsilonS(lshift);
  
  electron.Calculate(0);
  
  std::cout << std::setprecision(10);
  
  std::cout << "Optimized energy: " << electron.getEnergy() << " Hartree\n\n";

  //get info on FMOs
  std::vector<double> HOMOorb;
  double Ehomo = electron.getHOMO(HOMOorb);
  std::vector<double> LUMOorb;
  double Elumo = electron.getLUMO(LUMOorb);
  std::cout << "HOMO-LUMO gap " << Elumo - Ehomo << std::endl << std::endl;
  

//get charges
  //replace "Mulliken" with "Loewdin", "CM1", "CM2", "CM3", "CM5", "Goedecker"
  //careful becayse most CMx charges only available for AM1 and PM3
  auto charge_style = { "Mulliken",  "CM1", "CM2", "CM3", "CM5", "Goedecker"};
  for (auto style : charge_style)
  {
   std::cout << style << " Model\n";
   std::vector<double> charges = electron.getCharges(style);

   std::cout << style<<" charges:" << std::endl;
   for (size_t idAtm = 0; idAtm < charges.size(); ++idAtm) {
    std::cout << atoms[idAtm] << "    " << charges[idAtm] << std::endl;
   }
  }

  return 0;
}
