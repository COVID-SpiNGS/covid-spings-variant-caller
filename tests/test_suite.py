import unittest

from vc_queue_test import VCQueueTest

# SRC: https://stackoverflow.com/questions/31556718/is-it-possible-to-run-all-unit-test

def create_suite():
    test_suite = unittest.TestSuite()
    test_suite.addTest(VCQueueTest())
    return test_suite


if __name__ == '__main__':
    suite = create_suite()

    runner = unittest.TextTestRunner()
    runner.run(suite)
