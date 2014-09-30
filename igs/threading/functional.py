# Parallel implementations of various functions
import Queue

from igs.threading import threads
from time import sleep

def pmap(f, iterable, num_workers=1):
    def _worker(work_queue, result):
        while not work_queue.empty():
            idx, work = work_queue.get()
            sleep(0.1)
            result.append((idx, f(work)))
            sleep(0.1)
            work_queue.task_done()
        while True:
            try:
                if work_queue.empty() is True:
                    break
            except KeyboardInterrupt:
                work_queue = Queue.Queue()
                break
                
        
    # We want to ensure the order is the same
    # on the output string so we index each value
    # so we can reconstruct it
    work_queue = Queue.Queue()
    for idx, v in enumerate(iterable):
        try:
            work_queue.put((idx, v))
            sleep(0.1)
        except:
            pass

    results = []
    worker_threads = []
    for i in range(num_workers):
        results.append([])
        worker_threads.append(threads.runThread(_worker, work_queue, results[i]))
        sleep(0.1)

    result = []

    for r in results:
        result.extend(r)

    result.sort()

    for th in worker_threads:
        th.join()
    
    return [v for _, v in result]




    
