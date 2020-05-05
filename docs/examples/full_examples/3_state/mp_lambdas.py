import numpy as np
import biceps
import multiprocessing as mp

####### Data and Output Directories #######
energies = np.loadtxt('energies.dat')*627.509/0.5959 # from hartrees to reduced free energies F = f/kT
energies -= energies.min()  # set ground state to zero, just in case
states = len(energies)
input_data = biceps.toolbox.sort_data('NOE')
print(f"Input data: {biceps.toolbox.list_extensions(input_data)}")
####### Parameters #######
nsteps=1000000
n_lambdas = 2
nreplicas = 100
outdir = '%s_steps_%s_replica_%s_lam'%(nsteps, nreplicas, n_lambdas)
biceps.toolbox.mkdir(outdir)
print(f"nSteps of sampling: {nsteps}\nnReplicas: {nreplicas}")
lambda_values = np.linspace(0.0, 1.0, n_lambdas)
parameters = [dict(ref="uniform", uncern=(0.05, 5.0, 1.02), gamma=(0.2, 5.0, 1.02)),]
#parameters = [dict(ref="exp", uncern=(0.05, 5.0, 1.02), gamma=(0.2, 5.0, 1.02)),]
####### Multiprocessing Lambda values #######
def mp_lambdas(Lambda):
    ensemble = biceps.Ensemble(Lambda, energies)
    ensemble.initialize_restraints(input_data, parameters)
    sampler = biceps.PosteriorSampler(ensemble.to_list(), nreplicas,
            )
            #freq_write_traj=100., freq_save_traj=10.)
    sampler.sample(nsteps=nsteps, print_freq=1000, verbose=False)
    sampler.traj.process_results(outdir+'/traj_lambda%2.2f.npz'%(lam))
    filename = outdir+'/sampler_lambda%2.2f.pkl'%(lam)
    biceps.toolbox.save_object(sampler, filename)
    print('...Done.')
print("Number of CPU's: %s"%(mp.cpu_count()))
p = mp.Pool(processes=n_lambdas) # knows the number of CPU's to allocate
print(f"Number of processes: {n_lambdas}")
jobs = []
for lam in lambda_values:
    process = p.Process(target=mp_lambdas, args=(lam,))
    jobs.append(process)
    jobs[-1].start()
    active_processors = [jobs[i].is_alive() for i in range(len(jobs))]
    if (len(active_processors) == mp.cpu_count()-1) and all(active_processors) == True:
        while all(active_processors) == True:
            active_processors = [jobs[i].is_alive() for i in range(len(jobs))]
        inactive = int(np.where(np.array(active_processors) == False)[0])
        jobs[inactive].terminate()
        jobs.remove(jobs[inactive])
for job in jobs:
    job.join()
p.close()

####### Posterior Analysis #######
A = biceps.Analysis(states=states, resultdir=outdir,
    BSdir='BS.dat', popdir='populations.dat',
    picfile='BICePs.pdf')
A.plot()
