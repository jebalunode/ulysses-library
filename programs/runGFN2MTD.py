import os
import sys
import subprocess

#arguments and give them default values
exe = "./GFN-MTD.exe"
geometryfile = ""
charge = "0"
Telec = "300"
solvation = "0"
solventname = "gas"
optimisegeometry = "0"
SHAKE = "0"
equilibrate = "0"
maxtime = "500"          #ps
timestep = "0.001"       #ps
geomstorefreq = "0.01"   #ps
doMTD = "0"
MTDstruct = "1000"
kappa = "0.02"
alpha = "1.0"
MTDcollect = "100"
Temperature = "300"
driftthreshold = "2.0e-3"
printingfreq = "200"
mult = "1"
storefolder = "md/"
basename = "md"
mtdatoms = "mtd.txt"
biasgeom = ""

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
    elif argument.startswith("shake"):
        SHAKE = argument.replace("shake","")
    elif argument.startswith("equilibrate") or argument.startswith("do_equilibration") or argument.startswith("equil"):
        equilibrate = argument.replace("equilibrate","").replace("do_equilibration","").replace("equil","")
    elif argument.startswith("maxtime") or argument.startswith("maxsimtime") or argument.startswith("maximumtime"):
        maxtime = argument.replace("maxtime","").replace("maxsimtime","").replace("maximumtime","")
    elif argument.startswith("timestep") or argument.startswith("time_step") or argument.startswith("step"):
        timestep = argument.replace("timestep","").replace("time_step","").replace("step","")
    elif argument.startswith("storingfreq") or argument.startswith("geomstorefreq"):
        geomstorefreq = argument.replace("storingfreq","").replace("geomstorefreq","")
    elif argument.startswith("domtd") or argument.startswith("dometadynamics"):
        aux = argument.replace("domtd","").replace("dometadynamics","")
        if (aux == "true") or (aux == "1"): doMTD = "1"
    elif argument.startswith("mtdstruct") or argument.startswith("numbermtdstructures"):
        MTDstruct = argument.replace("mtdstruct","").replace("numbermtdstructures","")
    elif argument.startswith("kappa"):
        kappa = argument.replace("kappa","")
    elif argument.startswith("alpha"):
        alpha = argument.replace("alpha","")
    elif argument.startswith("mtdcollectionfreq") or argument.startswith("mtdcollect"):
        MTDcollect = argument.replace("mtdcollectionfreq","").replace("mtdcollect","")
    elif argument.startswith("temperature") or argument.startswith("temp"):
        Temperature = argument.replace("temperature","").replace("temp","")
    elif argument.startswith("drift"):
        driftthreshold = argument.replace("drift","")
    elif argument.startswith("printingfreq") or argument.startswith("printingfrequency"):
        printingfreq = argument.replace("printingfrequency","").replace("printingfreq","")
    elif argument.startswith("multiplicity") or argument.startswith("mult"):
        mult = argument.replace("multiplicity","").replace("mult","")
    elif argument.startswith("storefolder") or argument.startswith("folder2store"):
        storefolder = argument.replace("storefolder","").replace("folder2store","")
    elif argument.startswith("storefilename") or argument.startswith("basename"):
        basename = argument.replace("storefilename","").replace("basename","")
    elif argument.startswith("mtdrestraints") or argument.startswith("mtdatoms"):
        mtdatoms = argument.replace("mtdrestraints","").replace("mtdatoms","")
    elif argument.startswith("biasgeom") or argument.startswith("biasgeometry"):
        biasgeom = argument.replace("biasgeom","").replace("biasgeometry","")

#ensure consistency
if geometryfile + "blahblah" == "blahblah":
    print("ERROR: geometry file is needed")
    sys.exit()

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

if not storefolder.endswith("/"): storefolder += "/"

def RunMD():
    proc = subprocess.Popen([exe,geometryfile,charge,solvation,solventname,optimisegeometry,storefolder + basename,SHAKE,equilibrate,maxtime,timestep,geomstorefreq,doMTD,MTDstruct,kappa,alpha,MTDcollect,mtdatoms,Temperature,driftthreshold,printingfreq,mult,Telec,biasgeom],stdout = subprocess.PIPE,stderr = subprocess.PIPE)
    proc.wait()
    print(geometryfile,charge,solvation,solventname,optimisegeometry,storefolder + basename,SHAKE,equilibrate,maxtime,timestep,geomstorefreq,doMTD,MTDstruct,kappa,alpha,MTDcollect,mtdatoms,Temperature,driftthreshold,printingfreq,mult,Telec,biasgeom)
    out,err = proc.communicate()
    return out,err
out,err = RunMD()
