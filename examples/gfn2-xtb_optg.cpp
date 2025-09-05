#include "../src/GFN.hpp"
#include "../src/math/SolverPackage.hpp"
#include <stdlib.h>

//program to run GFN2-xTB calculations

int main(int argc, char** argv) {

  //argument 1 to program is the geometry in xyz format
  //argument 2 to program is the charge
  //argument 3 to program is the electronic temperature (default = 300 K)

  std::cout << "running " << argv[1] << "\n";

  char *p;
  int charge = strtol(argv[2],&p,10);

  char *q;
  double Telec = strtod(argv[3],&q);
  
  std::cout << "charge = " << charge << std::endl;
  std::cout << "T electron = " << Telec << std::endl;
  
  //allocation of molecule with certain charge and multiplicity 2S+1 = 1
  Molecule Mol1(argv[1],charge,1);
  
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
  
  //SCF, with no print (silent)
  electron.Calculate(0);
  
  //allocation of BFGS object, using the direct algorithm (works with the Hessian, not the inverted Hessian), using the dogleg
  //other line searches: replace number 6 with 
  //1 -> Davidon
  //2 -> Barzilai-Borwein
  //3 -> Armijo
  //4 -> More'-Thuente
  //5 -> Fletcher's parabollic interpolation
  //other types of Hessian update: replace the 4 with 
  //0 -> Murtagh-Sargent; B. A. Murtagh, R. W. H. Sargent, Comp. J., 13(2), 185, 1970
  //1 -> Powell-symmetric-Broyden
  //2 -> Bofill; J. M. Bofill, J. Comput. Chem., 15(1),1,1994
  //3 -> Bakken-Helgaker; V. Bakken, T. Helgaker, J. Chem. Phys., 117(20), 9160, 2002
  BFGSd solve(4,6);
  //other solvers: 
  //BakerRFO solve(4,0,0); //Baker's RFO method; 4 is a Hessian update, may take other values
  //BFGSi instead of BFGSd, exact same parameters
  //the number 4 is the type of Hessian used; currently you have the Lindh Hessian; 3 is for Schlegel and 0 is numerical
  double energy_threshold = 5.0e-6;
  double gradient_threshold = 1.0e-3;
  SolverOpt(electron,solve,4,0,energy_threshold,gradient_threshold);
  
  std::cout << std::setprecision(10);
  
  //get the molecule with optimized geometry
  Molecule Mol2 = electron.Component();
  
  std::vector<size_t> atoms = Mol2.Atoms();

  //write geometry to xyz
  std::string geomfilename = argv[1];
  size_t extpos = geomfilename.find(".xyz");
  size_t lengeom = geomfilename.size();
  geomfilename.replace(extpos,lengeom,"");
  Mol2.WriteXYZ(geomfilename,0);
  //this writes the optimized geometry to the file adenine0.xyz
  
  //other properties may come now here


  //get charges
  //replace "Mulliken" with "Loewdin", "CM1", "CM2", "CM3", "CM5", "Goedecker"
  //careful becayse most CMx charges only available for AM1 and PM3
  auto charge_style = { "Mulliken", "Loewdin", "CM1", "CM2", "CM3", "CM5", "Goedecker"};
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
