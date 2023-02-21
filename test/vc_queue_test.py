import unittest

import unittest

from client_server.vc_queue import VCQueue
from client_server.vc_exception import VCException


class VCQueueTest(unittest.TestCase):

    vc_q = VCQueue(3)
    def setUp(self) -> None:
        size = 3
        self.vc_q = VCQueue(size)

    def test_creation(self):
        with self.assertRaises(VCException):
            vc_q2 = VCQueue(11)

    def test_put1(self):
        self.vc_q.put('test')
        self.assertEqual(self.vc_q.length(), 1)

    def test_put2(self):
        size = 3
        vc_queue = VCQueue(size)
        vc_queue.put(('test1', 'test2'))
        self.assertEqual(vc_queue.length(), 1)

    def test_process(self):
        size = 3
        vc_queue = VCQueue(size)
        vc_queue.put(('test1', 'test2'))
        self.assertEqual(vc_queue.length(), 1)
        vc_queue.process()
        self.assertEqual(vc_queue.length(), 0)

    def test_length(self):
        size = 3
        vc_queue = VCQueue(size)
        self.assertEqual(vc_queue.length(), 0)
        vc_queue.put(('test1', 'test2'))
        self.assertEqual(vc_queue.length(), 1)

    def test_isempty1(self):
        size = 3
        vc_queue = VCQueue(size)
        self.assertTrue(vc_queue.is_empty())

    def test_isempty2(self):
        size = 3
        vc_queue = VCQueue(size)
        vc_queue.put(('test1', 'test2'))
        self.assertFalse(vc_queue.is_empty())


if __name__ == '__main__':
    unittest.main()
