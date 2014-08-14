# Parallel implementations of various functions
import Queue

from igs.threading import threads
from time import sleep

def pmap(f, iterable, num_workers=1):
    def _worker(work_queue, result):
        while not work_queue.empty():
            try:
                idx, work = work_queue.get()
                result.append((idx, f(work)))
                work_queue.task_done()
            except:
                break
        
    # We want to ensure the order is the same
    # on the output string so we index each value
    # so we can reconstruct it
    work_queue = Queue.Queue()
    for idx, v in enumerate(iterable):
        work_queue.put((idx, v))

    results = []
    worker_threads = []
    for i in range(num_workers):
        results.append([])
        worker_threads.append(threads.runThread(_worker, work_queue, results[i]))

    result = []
    while num_workers > 0:
        for r in results:
            result.extend(r)
            num_workers -= 1

    result.sort()

    for th in worker_threads:
        th.join()
    
    return [v for _, v in result]




    
