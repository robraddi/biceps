import os, sys, pickle
import numpy as np
import biceps
import multiprocessing as mp

####### Data and Output Directories #######
energies = np.loadtxt('cineromycin_B/cineromycinB_QMenergies.dat')*627.509  # convert from hartrees to kcal/mol
energies = energies/0.5959   # convert to reduced free energies F = f/kT
energies -= energies.min()  # set ground state to zero, just in case
states = len(energies)
print(f"Possible input data extensions: {biceps.toolbox.list_possible_extensions()}")
input_data = biceps.toolbox.sort_data('cineromycin_B/J_NOE')
print(f"Input data: {biceps.toolbox.list_extensions(input_data)}")
#outdir = '1M_steps_2_replica_2_lam'#'_results_test'
outdir = '1M_steps_6_replica_3_lam'

#outdir = '_results_test'
biceps.toolbox.mkdir(outdir)
####### Parameters #######
nsteps=1000000
maxtau = 1000
n_lambdas = 2
nreplicas = 4
lambda_values = np.linspace(0.0, 1.0, n_lambdas)
parameters = [
        dict(ref="uniform", uncern=(0.05, 20.0, 1.02)),
        dict(ref="exp", uncern=(0.05, 5.0, 1.02), gamma=(0.2, 5.0, 1.02)),
        ]
"""
for lam in lambda_values:
    print(f"lambda: {lam}")
    ensemble = biceps.Ensemble(lam, energies)
    ensemble.initialize_restraints(input_data, parameters)
    sampler = biceps.PosteriorSampler(ensemble.to_list(), nreplicas)
    sampler.sample(nsteps, print_freq=1, verbose=True)
    #sampler.sample(nsteps, verbose=True)
    sampler.traj.process_results(outdir+'/traj_lambda%2.2f.npz'%(lam))
    outfilename = 'sampler_lambda%2.2f.pkl'%(lam)
    fout = open(os.path.join(outdir, outfilename), 'wb')
    pickle.dump(sampler, fout)
    fout.close()
    print('...Done.')
    #exit()

"""

'''
####### Convergence Check #######
C = biceps.Convergence(trajfile=outdir+"/traj_lambda0.00.npz", resultdir=outdir)
C.get_autocorrelation_curves(method="normal", maxtau=maxtau)
C.plot_auto_curve(fname="auto_curve.pdf", xlim=(0, maxtau))
C.process(nblock=5, nfold=10, nround=100, savefile=True,
    plot=True, block=True, normalize=True)

#'''


####### Posterior Analysis #######
A = biceps.Analysis(states=states, resultdir=outdir,
    BSdir='BS.dat', popdir='populations.dat',
    picfile='BICePs.pdf')
A.plot()






