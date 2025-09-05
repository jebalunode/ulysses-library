#program to split a protein in a pdb file and run a series of calculations
import os
import sys
import subprocess

#arguments and give them default values
pdbfile = ""                    #pdb file to process
gfolder = "geometries/"         #folder storing the geometries of pairs
ofolder = "output/"             #folder storing the output of calculations
if not os.path.exists(gfolder): os.mkdir(gfolder)
if not os.path.exists(ofolder): os.mkdir(ofolder)
charge = 0                      #ligand charge
userdefcharg = False            #whether user defined charge
chain = "A"                     #chain identifier
cap = True                      #whether to cap the end of residues
ligandname = "0"                #define the name of the ligand
pH = 7.0                        #pH for protonation
listactiveres = []              #list of active residues, for which calculations are ran
activeres = "all"               #user supplied list
                                #this can be given in the form "1-45,47,48;60"
alpha = 1.0                     #alpha for in-pocket
kappa = 0.025                   #kappa for in-pocket
ignh = 1                        #whether to ignore protons
Telec = 300.0                   #electronic temperature
solvation = "1"                 #whether to use ALPB solvation
solventname = "water"           #solvent to use
calcdensity = "0"               #whether to calculate density
energythresh = "1.0e-6"         #threshold to use for energy convergence
gradientthresh = "1.0e-3"       #threshold to use for gradient convergence
do_inpocket = 1                 #whether to do in-pocket or uncontrained optimization

#parse whatever we get
for iarg in range(1,len(sys.argv)):
    argument = sys.argv[iarg].lower().replace("=","")
    if argument.startswith("pdb") or argument.startswith("pdbfile"):
        pdbfile = argument.replace("pdbfile","").replace("pdb","") + "pdb"
    elif argument.startswith("charge") or argument.startswith("chrg"):
        charge = int(argument.replace("charge","").replace("chrg",""))
        userdefcharg = True
    elif argument.startswith("chain"):
        chain = argument.replace("chain","")
    elif argument.startswith("telec") or argument.startswith("telectronic"):
        Telec = float(argument.replace("telectronic","").replace("telec",""))
    elif argument.startswith("cap"):
        aux = argument.replace("cap","")
        if (aux == "false") or (aux == "0"): cap = False
    elif argument.startswith("ignh") or argument.startswith("ignoreH") or argument.startswith("ignore_protons"):
        aux = argument.replace("ignh","").replace("ignoreH","").replace("ignore_protons","")
        if (aux == "false") or (aux == "0"): ignh = 0
    elif argument.startswith("doip") or argument.startswith("do_in-pocket") or argument.startswith("doinpocket") or argument.startswith("do_inpocket"):
        aux = argument.replace("doip","").replace("do_in-pocket","").replace("doinpocket","").replace("do_inpocket","")
        if (aux == "false") or (aux == "0"): do_inpocket = 0
    elif argument.startswith("alpha"):
        alpha = float(argument.replace("alpha",""))
    elif argument.startswith("kappa"):
        kappa = float(argument.replace("kappa",""))
    elif argument.startswith("ph"):
        pH = float(argument.replace("ph",""))
    elif argument.startswith("ligname") or argument.startswith("ligandname") or argument.startswith("ligand_name"):
        ligandname = argument.replace("ligname","").replace("ligandname","").replace("ligand_name","")
    elif argument.startswith("active_residues") or argument.startswith("activeresidues") or argument.startswith("activeres") or argument.startswith("ares"):
        activeres = argument.replace("active_residues","").replace("activeresidues","").replace("activeres","").replace("ares","")
    elif argument.startswith("usesolv") or argument.startswith("usesolvent"):
        solvation = argument.replace("usesolvent","").replace("usesolv","")
    elif argument.startswith("solvent") or argument.startswith("solventname"):
        solventname = argument.replace("solventname","").replace("solvent","")
    elif argument.startswith("calcdens") or argument.startswith("calcdensity") or argument.startswith("calculatedensity"):
        calcdensity = argument.replace("calcdensity","").replace("calcdens","").replace("calculatedensity","")
    elif argument.startswith("ethresh") or argument.startswith("energythresh") or argument.startswith("energythreshold"):
        energythresh = argument.replace("ethresh","").replace("energythreshold","").replace("energythresh","")
    elif argument.startswith("gthresh") or argument.startswith("gradthresh") or argument.startswith("gradientthreshold"):
        gradientthresh = argument.replace("gthresh","").replace("gradientthreshold","").replace("gradthresh","")

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

#ensure consistency, i.e., we have a pdb file
if pdbfile + "blahblah" == "blahblah":
    print("ERROR: geometry file is needed")
    sys.exit()

#functions
def Write2XYZ(atmlist, geometry, name):
    #function writing an xyz file given geometry and atom list
    wfile = open(name + ".xyz","w")
    wfile.write(str(len(atmlist)) + "\n\n")
    for idatm in range(len(atmlist)):
        wfile.write(atmlist[idatm] + "   " + str(geometry[idatm][0]) + "    " + str(geometry[idatm][1]) + "    " + str(geometry[idatm][2]) + "\n")
    wfile.close()

def CleanTheLine(line):
    #function that removes all double spaces from string
    cleanline = line
    for idx in range(250):
        if "  " in cleanline: 
            line = cleanline.replace("  "," ")
            cleanline = line
        else: break
    return cleanline

def Write2File(etotal, rmsd, totalpol, atmlist, charges, polarisabilities, converged, outfile):
    wfileres = open(outfile,"w")
    wfileres.write("Energy: " + str(etotal) + " Eh\n")
    wfileres.write("RMSD with starting structure: " + str(rmsd) + " Angstrom\n")
    wfileres.write("Optimisation Converged? " + str(converged) + "\n")
    wfileres.write("Atomic properties: charge    polarisability\n")
    totalcharge = 0.0
    for idatm in range(len(atmlist)):
        wfileres.write("    " + str(atmlist[idatm]) + "   " + str(charges[idatm]) + "   " + str(polarisabilities[idatm]) + "\n")
        totalcharge += charges[idatm]
    wfileres.write("Total Charge: " + str(totalcharge) + "\n")
    wfileres.write("Total Polarisability: " + str(totalpol) + " a0^3\n")
    wfileres.close()
    return 0

def ObabelProtonation(path2xyz):
    #function protonating the ligand using obabel
    #modify here to use pybel or for that matter of fact even RDkit
    #->
    intermedmol2 = path2xyz.replace(".xyz",".mol2")
    finalmol2 = intermedmol2.replace(".mol2","H.mol2")
    subprocess.run(["obabel",path2xyz,"-O",intermedmol2])
    subprocess.run(["obabel",intermedmol2,"-O",finalmol2,"-p",str(pH)])
    subprocess.run(["rm",intermedmol2])
    #<-
    #get the total charge of the molecule
    totalcharge = 0.0
    rfile = open(finalmol2,"r")
    rfilelines = rfile.readlines()
    fetch = False
    for line in rfilelines:
        if "@<TRIPOS>ATOM" in line: 
            fetch = True
            continue
        elif "@<TRIPOS>UNITY_ATOM_ATTR" in line: break
        elif "@<TRIPOS>BOND" in line: break
        if fetch:
            cleanline = CleanTheLine(line)
            data = cleanline.strip().rstrip().split(" ")
            totalcharge += float(data[8])
    return round(totalcharge)

def MOL2toXYZ(mol2file, restp):
    #conversion between mol2 and xyz
    #a few fixes are used for older versions of obabel that do not consider the pH
    geometry = []
    atoms = []
    rfile = open(mol2file,"r")
    rfilelines = rfile.readlines()
    fetch = False
    CDCZ = []
    NE = []
    for line in rfilelines:
        if line.startswith("@<TRIPOS>ATOM"):
            fetch = True
            continue
        elif "@<TRIPOS>UNITY_ATOM_ATTR" in line: break
        elif line.startswith("@<TRIPOS>BOND"): break
        if fetch:
            cleanline = CleanTheLine(line)
            data = cleanline.strip().rstrip().split(" ")
            aux = [float(data[2]),float(data[3]),float(data[4])]
            element = data[5]
            posdot = element.rfind(".")
            if posdot >= 0: element = data[5][:posdot]
            if restp.lower() == "glu" and ((data[1] == "HE1") or (data[1] == "HE2")): continue
            elif restp.lower() == "asp" and ((data[1] == "HD1") or (data[1] == "HD2")): continue
            elif restp.lower() == "arg" and ((data[1] == "CZ") or (data[1] == "CD")): CDCZ.append(aux)
            elif restp.lower() == "arg" and data[1] == "NE": NE = aux
            atoms.append(element)
            geometry.append(aux)
    #add missing protons
    if restp.lower() == "arg": 
        atoms.append("H")
        aux = [-0.5*(CDCZ[0][0] + CDCZ[1][0]),-0.5*(CDCZ[0][1] + CDCZ[1][1]),-0.5*(CDCZ[0][2] + CDCZ[1][2])]
        for idcoord in range(3):
            aux[idcoord] *= 1.5
            aux[idcoord] += 2.5*NE[idcoord]
        geometry.append(aux)
    Write2XYZ(atoms,geometry,mol2file.replace(".mol2",""))
    subprocess.run(["rm",mol2file])

def MergeXYZ(xyz1, xyz2, newname):
    #function that merges two geometry files into a new xyz file
    rfile = open(xyz1,"r")
    rfilelines = rfile.readlines()
    rfile.close()
    natm1 = int(rfilelines[0])
    atoms = []
    geometry = []
    for idatm in range(2,natm1 + 2):
        line = rfilelines[idatm]
        cleanline = CleanTheLine(line)
        data = cleanline.split(" ")
        atoms.append(data[0])
        geometry.append([float(data[1]),float(data[2]),float(data[3])])
    rfile = open(xyz2,"r")
    rfilelines = rfile.readlines()
    rfile.close()
    natm2 = int(rfilelines[0])
    for idatm in range(2,natm2 + 2):
        line = rfilelines[idatm]
        cleanline = CleanTheLine(line)
        data = cleanline.split(" ")
        atoms.append(data[0])
        geometry.append([float(data[1]),float(data[2]),float(data[3])])
    Write2XYZ(atoms,geometry,newname)

def RunIP(geometryfile,charge,alpha,kappa,ignH,ip,Telec,solvation,solventname,calcdensity,densfile,geomfile,energythresh,gradientthresh):
    #function to run the executable
    exefile = "./gfn2xtb_in_pocket_opt.exe"
    proc = subprocess.Popen([exefile,geometryfile,charge,alpha,kappa,ignH,ip,Telec,solvation,solventname,calcdensity,densfile,geomfile,energythresh,gradientthresh], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out,err = proc.communicate()
    return out,err

def ParseOutput(output):
    #function that parses the output of a calculation
    outputfile = output.decode("utf-8").split("\n")
    etotal = 0.0
    fetchatomicprops = False
    atmlist = []
    charges = []
    polarisabilities = []
    totalpol = 0.0
    rmsd = 0.0
    converged = True
    for iline in range(len(outputfile)):
        line = outputfile[iline]
        #check whether array data is available
        if line.startswith(">atom;charge;pol"): 
            fetchatomicprops = True
            continue
        elif line.startswith("<atom;charge;pol"): 
            fetchatomicprops = False
            continue
        #fetch data
        if line.startswith("Optimized g energy:"): etotal = float(line.replace("Optimized g energy:","").replace("Hartree","").replace(" ",""))
        elif line.startswith("RMSD with initial structure:"): rmsd = float(line.replace("RMSD with initial structure:","").replace(" ",""))
        elif "Total Polarizability" in line: totalpol = float(line.replace("Total Polarizability","").replace(" ",""))
        elif fetchatomicprops:
            cleanline = CleanTheLine(line)
            data = cleanline.split(" ")
            atmlist.append(int(data[0]))
            charges.append(float(data[2]))
            polarisabilities.append(float(data[3]))
        elif "no convergence in" in line: converged = False
    return etotal, rmsd, totalpol, atmlist, charges, polarisabilities, converged

#general variables of interest
geometry = []            #contains the geometry declared in the pdb files
atoms = []               #contains the list of atoms
atmtp = []               #list of atom types
residues = []            #contains the list of residues
resnumbr = []            #contains the residue numbering
chainid = []             #chain information
restp = []               #residue type, i.e., whether this is protein or something else

#process pdbfile
rfile = open(pdbfile,"r")
rfilelines = rfile.readlines()
for line in rfilelines:
    if line.startswith("ATOM") or line.startswith("HETATM"):
        atmnr = int(line[6:11])
        atmname = line[12:16].strip().rstrip()
        residue = line[16:20].strip().rstrip()
        if len(residue) == 4:
            #keep only one occupation, for now only A
            if not residue.startswith("A"): continue
            residue = residue[1:]
        schainid = line[21].strip().rstrip()
        resnr = int(line[22:26])
        xcoord = float(line[30:38])
        ycoord = float(line[38:46])
        zcoord = float(line[46:54])
        element = line[76:78].strip().rstrip()
        geometry.append([xcoord,ycoord,zcoord])
        atoms.append(element)
        residues.append(residue)
        atmtp.append(atmname)
        chainid.append(schainid)
        resnumbr.append(resnr)
        if line.startswith("ATOM"): restp.append("ATOM")
        elif line.startswith("HETATM"): restp.append("HETATM")
    if line.startswith("ENDMDL") or line.startswith("END"): break

#reduce to one chain, split protein and ligand
protein_geometry = []
protein_atoms = []
protein_atmtp = []
protein_residues = []
protein_resnumbr = []
ligand_geometry = []
ligand_atoms = []
for idatm in range(len(atoms)):
    if chainid[idatm] == chain:
        if restp[idatm] == "ATOM":
            #fetch for protein
            protein_geometry.append(geometry[idatm])
            protein_atoms.append(atoms[idatm])
            protein_atmtp.append(atmtp[idatm])
            protein_residues.append(residues[idatm])
            protein_resnumbr.append(resnumbr[idatm])
        elif restp[idatm] == "HETATM":
            #fetch for ligand
            if ligandname != "0": 
                if residues[idatm].lower() != ligandname: continue
            ligand_geometry.append(geometry[idatm])
            ligand_atoms.append(atoms[idatm])

#write base geometry of ligand
name4ligandfile = ligandname
if name4ligandfile == "0": name4ligandfile = "ligand"
Write2XYZ(ligand_atoms,ligand_geometry,gfolder + name4ligandfile)

#protonate the ligand and get its charge
obabelligcharge = ObabelProtonation(gfolder + name4ligandfile + ".xyz")
#write the final geometry to xyz
MOL2toXYZ(gfolder + name4ligandfile + "H.mol2","ligand")
subprocess.run(["rm",gfolder + name4ligandfile + ".xyz"])
if (obabelligcharge != charge): 
    if userdefcharg: 
        print("WARNING: obabel and user defined charge mismatch:",charge,obabelligcharge)
        print("                                    user defined:",charge)
        print("                                  obabel defined:",obabelligcharge)
        charge = int(input("please modify the file " + gfolder + name4ligandfile + "H.xyz to the correct protonation state and redefine the charge:"))
    else: charge = obabelligcharge

scharge = str(charge)    #charge as string

#convert active residues to list
if activeres == "all": 
    for ires in range(1,protein_resnumbr[len(protein_resnumbr) - 1] + 1): listactiveres.append(ires)
else:
    nares = activeres.replace(",",";")
    data = nares.split(";")
    for idt in range(len(data)):
        if "-" in data[idt]:
            limits = data[idt].split("-")
            llim = int(limits[0])
            ulim = int(limits[1]) + 1
            for ires in range(llim,ulim): listactiveres.append(ires)
        else: listactiveres.append(int(data[idt]))

#now store info on the active residues
aux_geom_ares = []
aux_atp_ares = []
aux_atoms_ares = []
aux_res_ares = []
#generate shape of arrays
for ires in range(protein_resnumbr[len(protein_resnumbr) - 1]):
    aux_geom_ares.append([])
    aux_atp_ares.append([])
    aux_atoms_ares.append([])
    aux_res_ares.append("0")
#populate
for iatm in range(len(protein_atoms)):
    ires = protein_resnumbr[iatm] - 1
    aux_geom_ares[ires].append(protein_geometry[iatm])
    aux_atp_ares[ires].append(protein_atmtp[iatm])
    aux_atoms_ares[ires].append(protein_atoms[iatm])
    aux_res_ares[ires] = protein_residues[iatm]
if cap:
    #if cap, then add the C=O on the N terminus, which is the first residue
    Nterminus = False
    for ires in range(len(aux_res_ares)):
        if aux_res_ares[ires] == "0": continue
        if not Nterminus: 
            Nterminus = True
            continue
        #look for the C and O atoms
        for idatm in range(len(aux_atp_ares[ires - 1])):
            if aux_atp_ares[ires - 1][idatm] == "C":
                aux_geom_ares[ires].append(aux_geom_ares[ires - 1][idatm])
                aux_atoms_ares[ires].append("C")
            elif aux_atp_ares[ires - 1][idatm] == "O":
                aux_geom_ares[ires].append(aux_geom_ares[ires - 1][idatm])
                aux_atoms_ares[ires].append("O")

#now clean
geom_ares = []
atp_ares = []
atoms_ares = []
res_ares = []
chr_ares = []
files_ares_residues = []
files_ares_residues_ligand = []
pair_files = []
for ires in listactiveres:
    geom_ares.append(aux_geom_ares[ires - 1])
    atp_ares.append(aux_atp_ares[ires - 1])
    atoms_ares.append(aux_atoms_ares[ires - 1])
    res_ares.append(aux_res_ares[ires - 1])
#and write
for ires in range(len(listactiveres)):
    if len(atoms_ares[ires]) == 0: continue
    Write2XYZ(atoms_ares[ires],geom_ares[ires],gfolder + res_ares[ires] + str(listactiveres[ires]))
    auxchr = ObabelProtonation(gfolder + res_ares[ires] + str(listactiveres[ires]) + ".xyz")
    #convert to xyz
    MOL2toXYZ(gfolder + res_ares[ires] + str(listactiveres[ires]) + "H.mol2",res_ares[ires])
    if res_ares[ires] == "ARG" and auxchr == 0: auxchr = 1
    chr_ares.append(auxchr)
    subprocess.run(["rm",gfolder + res_ares[ires] + str(listactiveres[ires]) + ".xyz"])
    #store the file
    files_ares_residues.append(gfolder + res_ares[ires] + str(listactiveres[ires]) + "H.xyz")
    #build the pair
    files_ares_residues_ligand.append(gfolder + ligandname + "_" + res_ares[ires] + str(listactiveres[ires]) + ".xyz")
    MergeXYZ(gfolder + name4ligandfile + "H.xyz",gfolder + res_ares[ires] + str(listactiveres[ires]) + "H.xyz",gfolder + ligandname + "_" + res_ares[ires] + str(listactiveres[ires]))

Eh2kcalmol = 627.5096080305927

#get the deformation energy of the ligand
out,err = RunIP(gfolder + name4ligandfile + "H.xyz",scharge,str(alpha),str(kappa),str(ignh),"1",str(Telec),str(solvation),solventname,str(calcdensity),gfolder + name4ligandfile + "_ipdens",gfolder + name4ligandfile + "_ip",str(energythresh),str(gradientthresh))
etotal_ip,rmsd_ip,totalpol_ip,atmlist,charges_ip,polarisabilities_ip,converged_ip = ParseOutput(out)
print("in-pocket on ligand done")
if not converged_ip: print("WARNING: in-pocket of ligand not converged")
out,err = RunIP(gfolder + name4ligandfile + "H.xyz",scharge,str(alpha),str(kappa),str(ignh),"0",str(Telec),str(solvation),solventname,str(calcdensity),gfolder + name4ligandfile + "_optdens",gfolder + name4ligandfile + "_opt",str(energythresh),str(gradientthresh))
etotal_opt,rmsd_opt,totalpol_opt,atmlist,charges_opt,polarisabilities_opt,converged_opt = ParseOutput(out)
print("unconstrained optimization on ligand done\n")
if not converged_opt: print("WARNING: unconstrained optimisation of ligand not converged")
print("Deformation Energy:",(etotal_ip - etotal_opt)*Eh2kcalmol,"kcal/mol")
print("RMSD in-pocket:",rmsd_ip,"Angstrom")
print("RMSD unconstrained:",rmsd_opt,"Angstrom")
Write2File(etotal_ip,rmsd_ip,totalpol_ip,atmlist,charges_ip,polarisabilities_ip,converged_ip,ofolder + "ligand_ip.out")
Write2File(etotal_opt,rmsd_opt,totalpol_opt,atmlist,charges_opt,polarisabilities_opt,converged_opt,ofolder + "ligand_opt.out")

#rename ligand file
os.rename(gfolder + name4ligandfile + "H.xyz",gfolder + name4ligandfile + ".xyz")

#now run in-pocket on residues and pairs
energy_residues = []
rmsd_residues = []
energy_residues_ligand = []
wfile = open("interaction_energies.csv","w")
wfile.write("residue;Eint(deformed,kcal/mol);RMSD(residue,A)\n")
for ires in range(len(listactiveres)):
    #start with residue
    out,err = RunIP(files_ares_residues[ires],str(chr_ares[ires]),str(alpha),str(kappa),str(ignh),str(do_inpocket),str(Telec),str(solvation),solventname,str(calcdensity),files_ares_residues[ires].replace(".xyz","_dens"),files_ares_residues[ires].replace(".xyz","_opt"),str(energythresh),str(gradientthresh))
    print("residue",ires + 1,"/",len(listactiveres),"done")
    #get data
    etotal_res,rmsd_res,totalpol_res,atmlist_res,charges_res,polarisabilities_res,converged_res = ParseOutput(out)
    #check and store
    if not converged_res: print("WARNING: optimisation of residue",res_ares[ires] + str(listactiveres[ires]),"not converged")
    Write2File(etotal_res,rmsd_res,totalpol_res,atmlist_res,charges_res,polarisabilities_res,converged_res,ofolder + res_ares[ires] + "_" + str(listactiveres[ires]) + ".out")
    #now the pair
    out,err = RunIP(files_ares_residues_ligand[ires],str(chr_ares[ires] + charge),str(alpha),str(kappa),str(ignh),str(do_inpocket),str(Telec),str(solvation),solventname,str(calcdensity),files_ares_residues_ligand[ires].replace(".xyz","_dens"),files_ares_residues_ligand[ires].replace(".xyz","_opt"),str(energythresh),str(gradientthresh))
    print("pair",ires + 1,"/",len(listactiveres),"done")
    etotal_pair,rmsd_pair,totalpol_pair,atmlist_pair,charges_pair,polarisabilities_pair,converged_pair = ParseOutput(out)
    #check and store
    if not converged_pair: print("WARNING: optimisation of ligand-residue",res_ares[ires] + str(listactiveres[ires]),"pair not converged")
    Write2File(etotal_pair,rmsd_pair,totalpol_pair,atmlist_pair,charges_pair,polarisabilities_pair,converged_pair,ofolder + res_ares[ires] + str(listactiveres[ires]) + "_lig_" + ".out")
    wfile.write(res_ares[ires] + str(listactiveres[ires]) + ";" + str((etotal_pair - etotal_ip - etotal_res)*Eh2kcalmol) + ";" + str(rmsd_res) + "\n")
wfile.close()
