#include <iostream>
#include "../src/MolecularDynamics.hpp"
#include "../src/GFN.hpp"

int main(int argc, char** argv) {

  //arguments
  // 0 exe
  // 1 - geometry file
  // 2 - charge
  // 3 - solvation?
  // 4 - solvent
  // 5 - optimise geometry?
  // 6 - files to store geometries along the MD
  // 7 - SHAKE: 0 no, 1 protons, 2 all
  // 8 - equilibrate?
  // 9 - maximum simulation time in ps
  //10 - time step in ps
  //11 - frequency (in ps) at which geometries are dumped
  //12 - do metadynamics?
  //13 - how many structures are kept in the metadynamics?
  //14 - metadynamics k, positive (repulsive) or negative (attractive)
  //15 - metadynamics alpha
  //16 - frequency with which metadynamics structures are collected
  //17 - file with MTD atom restraints
  //18 - temperature to run the MD
  //19 - threshold for drifting
  //20 - frequency for printing output
  //21 - multiplicity
  //22 - electronic temperature
  //23 - bias structure
  
  char *p;
  int charge = strtol(argv[2],&p,10);
  char *r;
  int solvation = strtol(argv[3],&r,10);
  char *v;
  int optgeom = strtol(argv[5],&v,10);
  bool optgeometry = (optgeom > 0);
  char *q;
  int SHAKE = strtol(argv[7],&q,10);
  char *s;
  int equilibrate = strtol(argv[8],&s,10);
  char *t;
  double tmax = strtod(argv[9],&t);          //500
  char *u;
  double tstep = strtod(argv[10],&u);        //0.001
  char *w;
  double toptg = strtod(argv[11],&w);        //0.01
  char *a;
  int dometadyn = strtol(argv[12],&a,10);
  char *b;
  int numbstruct = strtol(argv[13],&b,10);   //5
  char *c;
  double kappa = strtod(argv[14],&c);        //0.02
  char *d;
  double alpha = strtod(argv[15],&d);        //1.0
  char *e;
  int mtdcollect = strtol(argv[16],&e,10);   //10
  char *l;
  double rtemp = strtod(argv[18],&l);        //298
  char *m;
  double driftthresh = strtod(argv[19],&m);  //2e-3
  char *n;
  int printfreq = strtol(argv[20],&n,10);    //200
  char *y;
  int multiplicity = strtol(argv[21],&y,10);
  char *z;
  double Telec = strtod(argv[22],&z);
  
  //basics for the molecule
  Molecule Mol(argv[1],charge,multiplicity,"C1");
  BSet basis(Mol,"gfn2");
  GFN2 electron(basis,Mol);
  //use ALPB solvation?
  if (solvation > 0) {
    electron.setSolvent(argv[4]);
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
  
  //allocation of dynamics objects
  Dynamics MDobj(Mol,tmax,equilibrate,optgeometry,tstep,toptg);
  MDobj.setGeometryFile(argv[6]);
  MDobj.ApplyConstraints(false,SHAKE);
  MDobj.setIntegration("LeapFrog");
  
  //metadynamics?
  if (dometadyn > 0) {
    MDobj.setMetaDynamics(true,numbstruct,kappa,alpha,mtdcollect,0.03);
    std::vector<int> mtdrestr = getNumbersFromFile(argv[17]);
    for (size_t idatm = 0; idatm < mtdrestr.size(); ++idatm) {
      MDobj.addMetaAtom(mtdrestr[idatm]);
    }
    Molecule newmol(argv[23],charge,multiplicity,"C1");
    if ((newmol.Natoms() > 0)&&(kappa < 0.0)) {
      std::vector<matrixE> metaset;
      metaset.resize(numbstruct);
      matrixE geom = newmol.Geometry();
      for (size_t idmt = 0; idmt < numbstruct; ++idmt) {
        metaset[idmt] = geom;
      }
      MDobj.setMetaSet(metaset);
    }
  }
  
  MDobj.runMD(electron,rtemp,driftthresh,printfreq);
  
  return 0;
}
