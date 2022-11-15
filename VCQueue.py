from queue import Queue


class VCQueue:

    q = Queue(maxsize=5)

    def __init__(self, size: int):
        self.size = size
        self.q = Queue(maxsize=self.size)

    def put(self, action: str):
        self.q.put(action)

    def dequeue(self):
        rint("Queue is Full or Not:", q.full())
        print("Size of Queue:", q.qsize())
        print("Removing Elements:")
        print(q.get())
        print(q.get())
        print(q.get())
        print("Empty or Not??", q.empty())
        print(q.get())
        print("Empty or Not??", q.empty())
        print("Size of Queue:", q.qsize())

    def _write_vcf(path: str):
        logging.info(f'Writing VCF to {path}')
        liveVariantCaller.write_vcf(path)

    def process(self):
        if not self.q.empty():
            action = self.q.get()
            if action == 'process':

            elif action == 'write':

