from threading import Thread, Event
import psutil
import time
from datetime import datetime
import os
from pathlib import Path

LOG_FILENAME = 'memory.log'


class MonitorMemory(Thread):
    """
    Get the memory usage (RSS) of this python process + all the child processes, which should include all the calls to elastix etc.
    Write data to file.
    from: https://www.g-loaded.eu/2016/11/24/how-to-terminate-running-python-threads-using-signals/

    Usage
    -----
    mem = MonitorMemory(out_dir)
    # monitor memory until we stop it with mem.stop()

    """
    def __init__(self, out_dir: Path):
        super(MonitorMemory, self).__init__()
        self.log_file = out_dir / LOG_FILENAME
        self.shutdown_flag = Event()
        self.start()

    def run(self):

        with open(self.log_file, 'w') as fh:
            fh.write('time,mem %\n')

            while not self.shutdown_flag.is_set():

                current_process = psutil.Process(os.getpid())

                mem = current_process.memory_percent()

                for child in current_process.children(recursive=True):
                    mem += child.memory_percent()
                now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

                fh.write(f'{now}, {str(mem)}\n')

                time.sleep(2)

    def stop(self):
        self.shutdown_flag.set()
        self.join()

    def stopped(self):
        return self.shutdown_flag.is_set()
