import unittest

from client_server.live_server import VCServer


class VCServerTest(unittest.TestCase):
    def test_run(self):
        self.assertEqual(False, False)  # add assertion here
        vc_server = VCServer()


if __name__ == '__main__':
    unittest.main()
