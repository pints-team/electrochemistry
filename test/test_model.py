#
# Tests the electrochemistry models (have to be compiled first!)
#
import unittest
import math

DEFAULT = {
    'reversed': True,
    'Estart': 0.5,
    'Ereverse': -0.1,
    'omega': 9.0152,
    'phase': 0,
    'dE': 0.08,
    'v': -0.08941,
    't_0': 0.001,
    'T': 297.0,
    'a': 0.07,
    'c_inf': 1 * 1e-3 * 1e-3,
    'D': 7.2e-6,
    'Ru': 8.0,
    'Cdl': 20.0 * 1e-6,
    'E0': 0.214,
    'k0': 0.0101,
    'alpha': 0.53,
}

DEFAULT_POMS = {
    'reversed': False,
    'Estart': 0.6,
    'Ereverse': -0.1,
    #'omega': 60.05168,
    'omega': 6.05168,
    'phase': 0,
    'dE': 20e-3,
    'v': -0.1043081,
    't_0': 0.00,
    'T': 298.2,
    'a': 0.0707,
    'c_inf': 0.1 * 1e-3 * 1e-3,
    'Ru': 50.0,
    'Cdl': 0.000008,
    'Gamma': 0.7 * 53.0e-12,
    'alpha1': 0.5,
    'alpha2': 0.5,
    'E01': 0.368,
    'E02': 0.338,
    'E11': 0.227,
    'E12': 0.227,
    'E21': 0.011,
    'E22': -0.016,
    'k01': 7300,
    'k02': 7300,
    'k11': 1e4,
    'k12': 1e4,
    'k21': 2500,
    'k22': 2500,
    'alpha1': 0.5,
    'alpha2': 0.5,
    'alpha11': 0.5,
    'alpha12': 0.5,
    'alpha21': 0.5,
    'alpha22': 0.5
}


class TestModel(unittest.TestCase):

    def test_ec_unwrapped(self):
        """
        Runs a simple simulation
        """
        import electrochemistry
        import numpy as np
        # Create model
        model = electrochemistry.ECModel(DEFAULT)
        times = model.suggest_times()
        # Run simulation
        values = model.simulate(times)

        self.assertEqual(len(values), len(times))

    def test_poms_unwrapped(self):
        """
        Runs a simple simulation
        """
        import electrochemistry
        import numpy as np
        # Create model
        model = electrochemistry.POMModel(DEFAULT_POMS)
        times = model.suggest_times()
        values = model.simulate(times)

        self.assertEqual(len(values), len(times))



    def test_ec_wrapper(self):
        """
        Wraps a `pints.ForwardModel` around a model.
        """
        import electrochemistry
        import numpy as np

        # Create some toy data
        model = electrochemistry.ECModel(DEFAULT)
        times = model.suggest_times()
        values = model.simulate(times)

        # Test wrapper
        parameters = ['E0', 'k0', 'Cdl']
        pints_model = electrochemistry.PintsModelAdaptor(model, parameters)

        # Get real parameter values
        # Note: Retrieving them from ECModel to get non-dimensionalised form!
        real = np.array([model.params[x] for x in parameters])
        # Test simulation via wrapper class
        values2 = pints_model.simulate(real, times)
        self.assertEqual(len(values), len(values2))
        self.assertTrue(np.all(np.array(values) == np.array(values2)))


if __name__ == '__main__':
    unittest.main()
