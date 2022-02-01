# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import os # used for working with file paths
import csv
from Bio.Phylo.PAML import codeml

os.getcwd()
os.chdir('C:\\Users\\lexlu\\Desktop')

dataDir = "Genes"  # Alignments Directory
outDir = "Outputs" # Output Directory
tree = "4taxa.tree" # Tree File
outputfile = "PAML.csv" # Output sheet


outputfile = os.path.join(outDir, "PAML.csv")

with open(outputfile, 'wb') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerow(["Gene","Omega_for_each_branch","Downy","Hairy"])
    
    list_of_files = [s for s in os.listdir(dataDir) if s.endswith('.phy')]
    for filename in list_of_files:
        if filename[-4:]==".phy":
            print filename
            filepath = os.path.join(dataDir, filename)
            gene = filename.split(".")[0]
            alignment = filepath
            out_file1 = os.path.join(outDir,gene+"_model1.results")
            
            try:
                # model 1 (Null Model)

                cml = codeml.Codeml(alignment = alignment, tree = tree,
                    out_file = out_file1, working_dir = dataDir)

                cml.set_options(verbose=1)
                cml.set_options(CodonFreq=2)
                cml.set_options(cleandata=0)
                cml.set_options(fix_blength=0)
                cml.set_options(model=1)
                cml.set_options(NSsites=[0])
                cml.set_options(fix_omega=0)
                cml.set_options(omega=0)
                cml.set_options(clock=0)
                cml.set_options(ncatG=10)
                cml.set_options(runmode=0)
                cml.set_options(fix_kappa=0)
                cml.set_options(kappa=0.4)
                cml.set_options(fix_alpha=1)
                cml.set_options(alpha=0)
                cml.set_options(Malpha=0)
                cml.set_options(Small_Diff=1e-6)
                cml.set_options(method=1)
                cml.set_options(aaDist=0)
                cml.set_options(RateAncestor=0)
                cml.set_options(aaRatefile="wag.dat")
                cml.set_options(icode=0)
                cml.set_options(seqtype=1)
                cml.set_options(getSE=0)
                cml.set_options(noisy=9)
                cml.set_options(Mgene=0)
               
                run1 = cml.run(command="C:\\Program Files (x86)\\paml4.9j\\bin\\codeml",verbose = True)
    
                results = run1.get('NSsites')
                results_1 = results.get(0)
                w_branch = results_1.get('omega tree')
                w = results_1.get('parameters').get('omega')
                PP = w[4]
                PV = w[5]

                writer.writerow([gene,w_branch,PP,PV])

            except:
                pass
            
# For single file

outputfile = os.path.join(outDir, "all_loci_PAML.csv")
with open(outputfile, 'wb') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerow(["Gene","Omega_for_each_branch","Downy","Hairy"])
    
    filename= "all_loci.phy"
    filepath= "all_loci.phy"
    gene = filename.split(".")[0]
    alignment = filepath
    out_file1 = os.path.join(outDir,gene+"_model1.results")
    
    try:
        # model 1 (Null Model)

        cml = codeml.Codeml(alignment = alignment, tree = tree,
            out_file = out_file1, working_dir = dataDir)

        cml.set_options(verbose=1)
        cml.set_options(CodonFreq=2)
        cml.set_options(cleandata=0)
        cml.set_options(fix_blength=0)
        cml.set_options(model=1)
        cml.set_options(NSsites=[0])
        cml.set_options(fix_omega=0)
        cml.set_options(omega=0)
        cml.set_options(clock=0)
        cml.set_options(ncatG=10)
        cml.set_options(runmode=0)
        cml.set_options(fix_kappa=0)
        cml.set_options(kappa=0.4)
        cml.set_options(fix_alpha=1)
        cml.set_options(alpha=0)
        cml.set_options(Malpha=0)
        cml.set_options(Small_Diff=1e-6)
        cml.set_options(method=1)
        cml.set_options(aaDist=0)
        cml.set_options(RateAncestor=0)
        cml.set_options(aaRatefile="wag.dat")
        cml.set_options(icode=0)
        cml.set_options(seqtype=1)
        cml.set_options(getSE=0)
        cml.set_options(noisy=9)
        cml.set_options(Mgene=0)
       
        run1 = cml.run(command="C:\\Program Files (x86)\\paml4.9j\\bin\\codeml",verbose = True)

        results = run1.get('NSsites')
        results_1 = results.get(0)
        w_branch = results_1.get('omega tree')
        w = results_1.get('parameters').get('omega')
        PP = w[4]
        PV = w[5]

        writer.writerow([gene,w_branch,PP,PV])

    except:
        pass
