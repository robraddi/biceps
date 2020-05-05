import numpy as np
import biceps

####### Data and Output Directories #######
energies = np.loadtxt('energies.dat')*627.509/0.5959 # from hartrees to reduced free energies F = f/kT
energies -= energies.min()  # set ground state to zero, just in case
states = len(energies)
input_data = biceps.toolbox.sort_data('NOE')
print(f"Input data: {biceps.toolbox.list_extensions(input_data)}")
#outdir = 'reg_results'
outdir = 'reg_results_uniform'
biceps.toolbox.mkdir(outdir)
####### Parameters #######
nsteps=1000000
n_lambdas = 3
nreplicas = 1
print(f"nSteps of sampling: {nsteps}\nnReplicas: {nreplicas}")
lambda_values = np.linspace(0.0, 1.0, n_lambdas)
parameters = [dict(ref="exp", uncern=(0.05, 4.0, 1.02), gamma=(0.2, 4.0, 1.02)),]
#parameters = [dict(ref="uniform", uncern=(0.05, 4.0, 1.02), gamma=(0.2, 3.0, 1.02)),]
'''
for lam in lambda_values:
    print(f"lambda: {lam}")
    ensemble = biceps.Ensemble(lam, energies)
    ensemble.initialize_restraints(input_data, parameters)
    sampler = biceps.PosteriorSampler(ensemble.to_list(), nreplicas)
    #sampler.sample(nsteps, print_freq=1, verbose=True)
    sampler.sample(nsteps, verbose=True)
    sampler.traj.process_results(outdir+'/traj_lambda%2.2f.npz'%(lam))
    filename = outdir+'/sampler_lambda%2.2f.pkl'%(lam)
    biceps.toolbox.save_object(sampler, filename)
    print('...Done.')

'''
####### Posterior Analysis #######
A = biceps.Analysis(states=states, resultdir=outdir,
    BSdir='BS.dat', popdir='populations.dat',
    picfile='BICePs.pdf')
A.plot(show=True)
