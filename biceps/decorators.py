'''
Author: Rob Raddi
Edits on: May 13th, 2020
'''

import multiprocessing as mp
from multiprocessing import Manager
import numpy as np

def multiprocess(*args, **kwargs):
    """Decorator method for multiprocessing functions.
    Please refer to https://docs.python.org/3/library/multiprocessing.html#multiprocessing.sharedctypes.multiprocessing.Manager
    for manager types.

    https://docs.python.org/3.8/library/multiprocessing.html#contexts-and-start-methods
    # FIXME: Python >= 3.8
    #process = p.Process(target=function, args=(iter,), ctx=mp.get_context(method='fork'))
    """

    def wrapper(function):
        n = len(kwargs['iterable'])
        print("Number of CPUs: %s"%(mp.cpu_count()))
        p = mp.Pool(processes=mp.cpu_count())
        print(f"Number of processes: {n}")
        jobs = []
        for iter in kwargs['iterable']:
            process = p.Process(target=function, args=(iter,))
            # FIXME: Python >= 3.8
            #process = p.Process(target=function, args=(iter,), ctx=mp.get_context(method='fork'))
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
    return wrapper


# NOTE: An example of mine that will help with Hamiltonian replica exchange.
'''
with Manager() as manager:
    results = manager.list()
    cv = KFold(folds, True, 1)
    list_of_dict = [{"fold":f, "indices":i} for f,i in enumerate(cv.split(X))]
    @biceps.multiprocess(iterable=list_of_dict)
    def function(list_of_dict):
        fold, indices = list_of_dict["fold"], list_of_dict["indices"]
        results.append([fold,indices])

'''







