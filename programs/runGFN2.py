import os
import sys
import subprocess

#arguments and give them default values
exe = "./GFN2program.exe"
geometryfile = ""
charge = "0"
Telec = "300"
solvation = "0"
solventname = "gas"
optimisegeometry = "0"
newgeometryfile = ""
thermo = "0"
energythresh = "1.0e-6"
gradientthresh = "1.0e-3"
calcdensity = "0"
densityfile = ""
electronicreactivity = "0"
orbitalreactivity = "0"
KoopmanIP = "0"
IP = "0"
EA = "0"
electronegativity = "0"
hardness = "0"

#parse whatever we get
for iarg in range(1,len(sys.argv)):
    argument = sys.argv[iarg].lower().replace("=","")
    if argument.startswith("geometry") or argument.startswith("geom"):
        geometryfile = argument.replace("geometry","").replace("geom","")
    elif argument.startswith("charge") or argument.startswith("chrg"):
        charge = argument.replace("charge","").replace("chrg","")
    elif argument.startswith("telec") or argument.startswith("telectronic"):
        Telec = argument.replace("telectronic","").replace("telec","")
    elif argument.startswith("usesolv") or argument.startswith("usesolvent"):
        solvation = argument.replace("usesolvent","").replace("usesolv","")
    elif argument.startswith("solvent") or argument.startswith("solventname"):
        solventname = argument.replace("solventname","").replace("solvent","")
    elif argument.startswith("optg") or argument.startswith("optimisegeom") or argument.startswith("optimizegeom") or argument.startswith("optimisegeometry") or argument.startswith("optimizegeometry"):
        optimisegeometry = argument.replace("optg","").replace("optimise","").replace("optimize","").replace("geometry","").replace("geom","")
    elif argument.startswith("newgeom") or argument.startswith("newgeomfile") or argument.startswith("newgeometryfile"):
        newgeometryfile = argument.replace("newgeometryfile","").replace("newgeomfile","").replace("newgeom","")
    elif argument.startswith("dothermo") or argument.startswith("thermo") or argument.startswith("thermodynamics"):
        thermo = argument.replace("dothermo","").replace("thermodynamics","").replace("thermo","")
    elif argument.startswith("ethresh") or argument.startswith("energythresh") or argument.startswith("energythreshold"):
        energythresh = argument.replace("ethresh","").replace("energythreshold","").replace("energythresh","")
    elif argument.startswith("gthresh") or argument.startswith("gradthresh") or argument.startswith("gradientthreshold"):
        gradientthresh = argument.replace("gthresh","").replace("gradientthreshold","").replace("gradthresh","")
    elif argument.startswith("calcdens") or argument.startswith("calcdensity") or argument.startswith("calculatedensity"):
        calcdensity = argument.replace("calcdensity","").replace("calcdens","").replace("calculatedensity","")
    elif argument.startswith("densityfile") or argument.startswith("densfile"):
        densityfile = argument.replace("densityfile","").replace("densfile","")
    elif argument.startswith("ereact") or argument.startswith("electronicreactivity"):
        electronicreactivity = argument.replace("ereact","").replace("electronicreactivity","")
    elif argument.startswith("oreact") or argument.startswith("orbreact") or argument.startswith("orbitalreactivity"):
        orbitalreactivity = argument.replace("oreact","").replace("orbreact","").replace("orbitalreactivity","")
    elif argument.startswith("koopman") or argument.startswith("ipkoopman") or argument.startswith("koopmanip"):
        KoopmanIP = argument.replace("ip","").replace("koopman","")
    elif argument.startswith("ip") or argument.startswith("ionizationpotential") or argument.startswith("ionisationpotential") or argument.startswith("ipot"):
        IP = argument.replace("ip","").replace("ionizationpotential","").replace("ionisationpotential","").replace("ipot","")
    elif argument.startswith("ea") or argument.startswith("electronaffinity") or argument.startswith("eaffin"):
        EA = argument.replace("ea","").replace("electronaffinity","").replace("eaffin","")
    elif argument.startswith("electronegativity"):
        electronegativity = argument.replace("electronegativity","")
    elif argument.startswith("hardness"):
        hardness = argument.replace("hardness","")

#ensure consistency
if geometryfile + "blahblah" == "blahblah":
    print("ERROR: geometry file is needed")
    sys.exit()

newgeometryfile.replace(".xyz","").replace(".sdf","").replace(".mol2","").replace(".pdb","")
if (newgeometryfile + "blahblah" == "blahblah") and (optimisegeometry == "1"): newgeometryfile = geometryfile.replace(".xyz","").replace(".sdf","").replace(".mol2","").replace(".pdb","") + "_opt"

if (solvation == "1"):
    if (solventname == "water") or (solventname == "h2o") or (solventname == "o"): solventname = "water"
    elif (solventname == "acetone") or (solventname == "cc(o)c"): solventname = "acetone"
    elif (solventname == "acetonitrile") or (solventname == "ch3cn") or (solventname == "ccn"): solventname = "acetonitrile"
    elif (solventname == "aniline") or (solventname == "phnh2") or (solventname == "nc1ccccc1") or (solventname == "c1ccc(cc1)n"): solventname = "aniline"
    elif (solventname == "benzaldehyde") or (solventname == "phcho") or (solventname == "occ1ccccc1") or (solventname == "c1ccc(cc1)co"): solventname = "benzaldehyde"
    elif (solventname == "benzene") or (solventname == "c6h6") or (solventname == "phh") or (solventname == "c1ccccc1"): solventname = "benzene"
    elif (solventname == "dichloromethane") or (solventname == "ch2cl2") or (solventname == "c(cl)cl") or (solventname == "c(cl)(cl)"): solventname = "dichloromethane"
    elif (solventname == "chloroform") or (solventname == "chcl3") or (solventname == "c(cl)(cl)cl") or (solventname == "c(cl)(cl)(cl)"): solventname = "chloroform"
    elif (solventname == "carbon disulfide") or (solventname == "carbondisulfide") or (solventname == "cs2") or (solventname == "scs"): solventname = "carbon disulfide"
    elif (solventname == "dioxane") or (solventname == "o1ccocc1"): solventname = "dioxane"
    elif (solventname == "dmf") or (solventname == "dimethylformamide") or (solventname == "cn(c)co"): solventname = "dmf"
    elif (solventname == "dmso") or (solventname == "dimethylsulfoxide") or (solventname == "cs(o)c") or (solventname == "cs(c)o"): solventname = "dmso"
    elif (solventname == "ethanol") or (solventname == "etoh") or (solventname == "ch3ch2oh") or (solventname == "cco"): solventname = "ethanol"
    elif (solventname == "diethyl ether") or (solventname == "etoet") or (solventname == "ccocc") or (solventname == "ch3ch2och2ch3"): solventname = "diethyl ether"
    elif (solventname == "ethyl acetate") or (solventname == "acoet") or (solventname == "etoac") or (solventname == "ccoc(o)c"): solventname = "ethyl acetate"
    elif (solventname == "furan") or (solventname == "furane"): solventname = "furane"
    elif (solventname == "hexadecane") or (solventname == "c16"): solventname = "hexadecane"
    elif (solventname == "hexane") or (solventname == "cccccc") or (solventname == "c6"): solventname = "hexane"
    elif (solventname == "methanol") or (solventname == "meoh") or (solventname == "co") or (solventname == "ch3oh"): solventname = "methanol"
    elif (solventname == "nitromethane") or (solventname == "meno2") or (solventname == "ch3no2") or (solventname == "cn(o)o"): solventname = "nitromethane"
    elif (solventname == "octanol") or (solventname == "cccccccc") or (solventname == "c8"): solventname = "octanol"
    elif (solventname == "phenol") or (solventname == "phoh") or (solventname == "oc1ccccc1") or (solventname == "c1ccc(cc1)o"): solventname = "phenol"
    elif (solventname == "thf") or (solventname == "tetrahydrofuran"): solventname = "thf"
    elif (solventname == "toluene") or (solventname == "phme") or (solventname == "cc1ccccc1"): solventname = "toluene"
    elif (solventname == "octanol wet") or (solventname == "octanolwet") or (solventname == "wet octanol") or (solventname == "wetoctanol"): solventname = "octanol wet"

if (densityfile + "blahblah" == "blahblah") and (calcdensity == "1"): densityfile = geometryfile.replace(".xyz","").replace(".sdf","").replace(".mol2","").replace(".pdb","") + "_dens"

def RunQMOptg():
    proc = subprocess.Popen([exe,geometryfile,charge,Telec,solvation,solventname,optimisegeometry,newgeometryfile,thermo,energythresh,gradientthresh,calcdensity,densityfile,electronicreactivity,orbitalreactivity,KoopmanIP,IP,EA,electronegativity,hardness], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out,err = proc.communicate()
    return out,err
out,err = RunQMOptg()

def ReadXYZ(filename):
    #function reading xyz format
    rfile = open(filename,"r")
    rfilelines = rfile.readlines()
    natoms = int(rfilelines[0])
    geometry = []
    for iatm in range(natoms):
        line = rfilelines[iatm + 2]
        cleanline = line
        for idx in range(250):
            line = cleanline.replace("  "," ")
            cleanline = line
            if not ("  " in cleanline): break
        data = cleanline.split(" ")
        aux = [float(data[1]),float(data[2]),float(data[3])]
        geometry.append(aux)
    return geometry

if err.decode("utf-8") + "blahblah" != "blahblah": print("error:",err)
else:
    #process data from output file
    outputfile = out.decode("utf-8").split("\n")
    print(out.decode("utf-8"))
    etotal = 0.0
    gsolvation = 0.0
    fetchvibr = False
    vibrations = []
    fetchthermo = False
    temperature = []
    entropy = []
    enthalpy = []
    gibbs = []
    internalenergy = []
    helmholtz = []
    cp = []
    cv = []
    fetchatomicprops = False
    atmlist = []
    charges = []
    polarisabilities = []
    totalpol = 0.0
    fetchelecreact = False
    fetchorbreact = False
    eFukuiElctrophilie = []
    eFukuiNucleophilie = []
    eFukuiRadical = []
    eSoftnessElctrophilie = []
    eSoftnessNucleophilie = []
    eSoftnessRadical = []
    oFukuiElctrophilie = []
    oFukuiNucleophilie = []
    oFukuiRadical = []
    oSoftnessElctrophilie = []
    oSoftnessNucleophilie = []
    oSoftnessRadical = []
    ipkoopman = 0.0
    geometry = []
    ip = 0.0
    ea = 0.0
    eneg = 0.0
    hard = 0.0
    for iline in range(len(outputfile)):
        line = outputfile[iline]
        #check whether array data is available
        if line.startswith(">all vibrational frequencies"): 
            fetchvibr = True
            continue
        elif line.startswith("<all vibrational frequencies"): 
            fetchvibr = False
            continue
        elif line.startswith(">Thermodynamics"): 
            fetchthermo = True
            continue
        elif line.startswith("<Thermodynamics"): 
            fetchthermo = False
            continue
        elif line.startswith(">atom;charge;pol"): 
            fetchatomicprops = True
            continue
        elif line.startswith("<atom;charge;pol"): 
            fetchatomicprops = False
            continue
        elif line.startswith(">Electronic Reactivity indices"): 
            fetchelecreact = True
            continue
        elif line.startswith("<Electronic Reactivity indices"): 
            fetchelecreact = False
            continue
        elif line.startswith(">Orbital Reactivity indices"): 
            fetchorbreact = True
            continue
        elif line.startswith("<Orbital Reactivity indices"): 
            fetchorbreact = False
            continue
        #fetch data
        if line.startswith("Total Energy = "): etotal = float(line.replace("Total Energy = ",""))
        elif line.startswith("Gsolv = "): gsolvation = float(line.replace("Gsolv = ",""))
        elif "Total Polarizability" in line: totalpol = float(line.replace("Total Polarizability","").replace(" ",""))
        elif line.startswith("Ionization Potential (Koopman):"): ipkoopman = float(line.replace("Ionization Potential (Koopman):","").replace("   eV",""))
        elif line.startswith("Ionization Potential (Definition): "): ip = float(line.replace("Ionization Potential (Definition): ","").replace("   eV",""))
        elif line.startswith("Electron Affinity (Definition): "): ea = float(line.replace("Electron Affinity (Definition): ","").replace("   eV",""))
        elif line.startswith("Electronegativity: "): eneg = float(line.replace("Electronegativity: ","").replace("   eV",""))
        elif line.startswith("Hardness: "): hard = float(line.replace("Hardness: ","").replace("   eV",""))
        elif fetchvibr: vibrations.append(float(line))
        elif fetchthermo: 
            data = line.split(";")
            temperature.append(float(data[0]))
            entropy.append(float(data[1]))
            enthalpy.append(float(data[2]))
            gibbs.append(float(data[3]))
            internalenergy.append(float(data[4]))
            helmholtz.append(float(data[5]))
            cp.append(float(data[6]))
            cv.append(float(data[7]))
        elif fetchatomicprops:
            data = line.split(";")
            atmlist.append(int(data[0]))
            charges.append(float(data[1]))
            polarisabilities.append(float(data[2]))
        elif fetchelecreact:
            data = line.strip().rstrip().replace("  "," ").split(" ")
            eFukuiElctrophilie.append(float(data[0]))
            eFukuiNucleophilie.append(float(data[1]))
            eFukuiRadical.append(float(data[2]))
            eSoftnessElctrophilie.append(float(data[3]))
            eSoftnessNucleophilie.append(float(data[4]))
            eSoftnessRadical.append(float(data[5]))
        elif fetchorbreact:
            data = line.strip().rstrip().replace("  "," ").split(" ")
            oFukuiElctrophilie.append(float(data[0]))
            oFukuiNucleophilie.append(float(data[1]))
            oFukuiRadical.append(float(data[2]))
            oSoftnessElctrophilie.append(float(data[3]))
            oSoftnessNucleophilie.append(float(data[4]))
            oSoftnessRadical.append(float(data[5]))
    if (optimisegeometry == "0"): newgeometryfile = geometryfile
    if not newgeometryfile.endswith(".xyz"): newgeometryfile += ".xyz"
    geometry = ReadXYZ(newgeometryfile)
    #data collected; if there is a whether, then variable is boolean
    #geometryfile               the geometry of the molecule
    #charge                     the total charge
    #solventname                the name of the solvent
    #optimisegeometry           whether geometry was optimised
    #newgeometryfile            the name of the file containing optimised geometry
    #geometry                   contains the best geometry (either the input one if not optg takes place, or the optimized geometry if optg took place)
    #thermo                     whether thermodynamics was done
    #calcdensity                whether the electronic density was calculated
    #densityfile                the name of the file containing electronic density
    #electronicreactivity       whether electronic reactivity data was calculated
    #orbitalreactivity          whether orbital reacticity data was calculated
    #KoopmanIP                  whether the Koopman ionisation potential was calculated
    #IP                         whether the ionisation potential was calculated
    #EA                         whether electron affinity was calculated
    #electronegativity          whether electronegativity was calculated
    #hardness                   whether hardness was calculated
    #etotal                     the molecule's total energy in Hartree
    #gsolvation                 solvation Gibbs free energy in Hartree
    #vibrations                 array containing vibrational frequencies in cm-1
    #temperature                array containing temperatures
    #entropy                    array containing entropies
    #enthalpy                   array containing enthalpies
    #gibbs                      array containing Gibbs free energies
    #internalenergy             array containing internal energies
    #helmholtz                  array containing Helmholts energies
    #cp                         array containing heat capacities at constant pressure
    #cv                         array containing heat capacities at constant volume
    #atmlist                    list of atoms
    #charges                    partial charges in electrons
    #polarisabilities           atomic polarisabilities in bohr
    #totalpol                   total polarisability in cubic bohr
    #eFukuiElctrophilie         Fukui index for electrophilicity (calculated from charges)
    #eFukuiNucleophilie         Fukui index for nucleophilicity (calculated from charges)
    #eFukuiRadical              Fukui index for radical reactivity (calculated from charges)
    #eSoftnessElctrophilie      softness index for electrophilicity (calculated from charges)
    #eSoftnessNucleophilie      softness index for nucleophilicity (calculated from charges)
    #eSoftnessRadical           softness index for radical reactivity (calculated from charges)
    #oFukuiElctrophilie         Fukui index for electrophilicity (calculated from frontier orbitals)
    #oFukuiNucleophilie         Fukui index for nucleophilicity (calculated from frontier orbitals)
    #oFukuiRadical              Fukui index for radical reactivity (calculated from frontier orbitals)
    #oSoftnessElctrophilie      softness index for electrophilicity (calculated from frontier orbitals)
    #oSoftnessNucleophilie      softness index for nucleophilicity (calculated from frontier orbitals)
    #oSoftnessRadical           softness index for radical reactivity (calculated from frontier orbitals)
    #ipkoopman                  Koopman ionisation potential in eV
    #ip                         ionisation potential in eV
    #ea                         electron affinity in eV
    #eneg                       electronegativity in eV
    #hard                       hardness in eV
