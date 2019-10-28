#!/usr/bin/python
#Analytic Fast Ion Distribution function:
# Automate fitting of NUBEAM data to generate splined differentiable
# distribution functions of the particle constants of motion.

from sys import argv
import os

if len(argv) > 1:
    jobid = argv[1]

    # Condense multi-CPU ptcl data in [jobid]_debug_nbi_ptcl_state_cpui_n.cdf
    #  and state data in [jobid]_ps_ts1_state.cdf into a single cdf file.
    # Requires condense utility in exec path.
    cfilename = jobid+"_condensed.cdf"
    if not os.access(cfilename, os.R_OK):
        os.system("condense "+jobid)
        os.system("particle_sort "+cfilename)

    # Construct Jacobian by extensively sampling NUBEAM distribution.
    # Needs state data and condensed output from previous command.
    # Requires mcgen in exec path.
    jfilename = jobid+"_jacobian.cdf"
    if not os.access(jfilename, os.R_OK):
        # First create input file with namelists specifying file root, n_particles
        npjac = 16000000
        with open("mcgen.in", "w") as file_object:
            file_object.write("&strings froot='"+jobid+"' /\n")
            file_object.write("&params nparts="+str(npjac)+" /\n")
        os.system("mcgen")
        os.system("particle_sort "+jfilename)

    # Fit splines to the generated distributions.
    # Needs condensed particle data and jacobian from two prior commands.
    # Requres fitjac in exec path.
    if not (os.access("pdist-1.spl", os.R_OK) and
            os.access("pdist01.spl", os.R_OK)):
        os.system("fitjac "+jobid+" -log")

else: # No argument provided -> print usage message.
    print("Usage: "+argv[0]+" [NUBEAM job id]\n")

