#
# Tests if the electrochemistry stuff can be loaded without issues.
#
import unittest


class TestBasics(unittest.TestCase):

    def test_import(self):
        import electrochemistry

if __name__ == '__main__':
    unittest.main()
