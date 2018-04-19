from __future__ import print_function
from math import sqrt, pi
from electrochemistry import seq_electron_transfer3_explicit, e_implicit_exponential_mesh
import pints
import numpy as np
import copy


class ECModel:

    """Represents one electron transfer model in solution
            A + e- <-> B

    Args:
        params (dict): dictionary of parameters, containing these keys
                 'reversed', 'Estart','Ereverse','omega','phase','dE','v','T','a','c_inf','D'
                  'Ru',
                 'Cdl',
                 'E0',
                 'k0',
                 'alpha'
    """

    def __init__(self, params):
        try:
            print('creating ECModel with (dimensional) parameters:')
            print('\treversed: ', params['reversed'])
            print('\tEstart: ', params['Estart'])
            print('\tEreverse: ', params['Ereverse'])
            print('\tomega: ', params['omega'])
            print('\tphase: ', params['phase'])
            print('\tamplitude: ', params['dE'])
            print('\tdc scan rate: ', params['v'])
            print('\ttemperature: ', params['T'])
            print('\telectrode area: ', params['a'])
            print('\tc_inf: ', params['c_inf'])
            print('\tdiffusion constant: ', params['D'])
            print('\tRu: ', params['Ru'])
            print('\tCdl: ', params['Cdl'])
            print('\tE0: ', params['E0'])
            print('\tk0: ', params['k0'])
            print('\talpha: ', params['alpha'])
        except NameError as e:
            print('NameError: ', e.value)

        self.dim_params = copy.copy(params)

        if (params['reversed']):
            self.dim_params['E0'] = self.dim_params['Estart'] - (
                self.dim_params['E0'] - self.dim_params['Ereverse'])
            self.dim_params['Estart'], self.dim_params[
                'Ereverse'] = self.dim_params['Ereverse'], self.dim_params['Estart']
            self.dim_params['v'] = -self.dim_params['v']

        self.E0, self.T0, self.L0, self.I0 = self._calculate_characteristic_values(
        )

        self.params = dict()
        self.params['Estart'] = self.dim_params['Estart'] / self.E0
        self.params['Ereverse'] = self.dim_params['Ereverse'] / self.E0
        self.params['omega'] = 2 * pi * self.dim_params['omega'] * self.T0
        if self.dim_params['reversed']:
            self.params['phase'] = self.dim_params['phase'] + pi
        else:
            self.params['phase'] = self.dim_params['phase']
        self.params['dE'] = self.dim_params['dE'] / self.E0

        self.params['k0'] = self.dim_params[
            'k0'] * self.L0 / self.dim_params['D']
        self.params['alpha'] = self.dim_params['alpha']
        self.params['E0'] = self.dim_params['E0'] / self.E0
        self.params['Ru'] = self.dim_params['Ru'] * abs(self.I0) / self.E0
        self.params['Cdl'] = self.dim_params['Cdl'] * \
            self.dim_params['a'] * self.E0 / (abs(self.I0) * self.T0)

        self._nondim_params = {}
        self._nondim_params['Estart'] = self.params['Estart']
        self._nondim_params['Ereverse'] = self.params['Ereverse']
        self._nondim_params['omega'] = self.params['omega']
        self._nondim_params['phase'] = self.params['phase']
        self._nondim_params['dE'] = self.params['dE']
        self._nondim_params['k0'] = self.params['k0']
        self._nondim_params['alpha'] = self.params['alpha']
        self._nondim_params['E0'] = self.params['E0']
        self._nondim_params['Ru'] = self.params['Ru']
        self._nondim_params['Cdl'] = self.params['Cdl']

    def dimensionalise(self, value, name):
        if name == 'Estart':
            return value * self.E0
        elif name == 'Ereverse':
            return value * self.E0
        elif name == 'omega':
            return value / (2 * pi * self.T0)
        elif name == 'phase':
            return value
        elif name == 'dE':
            return value * self.E0
        elif name == 'k0':
            return value * self.dim_params['D'] / self.L0
        elif name == 'alpha':
            return value
        elif name == 'E0':
            return value * self.E0
        elif name == 'Ru':
            return value * self.E0 / self.I0
        elif name == 'Cdl':
            return value * self.I0 * self.T0 / (self.dim_params['a'] * self.E0)
        else:
            return NaN

    def non_dimensionalise(self, value, name):
        if name == 'Estart':
            return value / self.E0
        elif name == 'Ereverse':
            return value / self.E0
        elif name == 'omega':
            return value * (2 * pi * self.T0)
        elif name == 'phase':
            return value
        elif name == 'dE':
            return value / self.E0
        elif name == 'k0':
            return value / self.dim_params['D'] * self.L0
        elif name == 'alpha':
            return value
        elif name == 'E0':
            return value / self.E0
        elif name == 'Ru':
            return value / self.E0 * self.I0
        elif name == 'Cdl':
            return value / self.I0 / self.T0 * (self.dim_params['a'] * self.E0)
        else:
            return NaN

    def suggest_times(self):
        final_time = self.params['Ereverse'] - self.params['Estart']
        return np.linspace(0, final_time, 1000, dtype='double')

    def simulate(self, times):
        times = np.asarray(times, dtype='double')
        current = np.empty_like(times)
        e_implicit_exponential_mesh(params, current, times)
        return current, times

    def set_params_from_vector(self, vector, names):
        for value, name in zip(vector, names):
            self.params[name] = value
            self.dim_params[name] = self.dimensionalise(value, name)

    def get_params_from_vector(self, names):
        vector = np.zeros(len(names))
        for i in range(len(names)):
            vector[i] = self.params[names[i]]
        return vector

    def _calculate_characteristic_values(self):

        v = self.dim_params['v']
        T = self.dim_params['T']
        a = self.dim_params['a']

        # Faraday constant (C mol-1)
        F = 96485.3328959
        # gas constant (J K-1 mol-1)
        R = 8.314459848

        E_0 = R * T / F
        T_0 = abs(E_0 / v)

        D = self.dim_params['D']
        L_0 = sqrt(D * T_0)
        c_inf = self.dim_params['c_inf']

        if self.dim_params['reversed']:
            I_0 = -D * F * a * c_inf / L_0
        else:
            I_0 = D * F * a * c_inf / L_0

        return E_0, T_0, L_0, I_0


class POMModel:

    def __init__(self, params):
        try:
            print('creating POMModel with (dimensional) parameters:')
            print('\tEstart: ', params['Estart'])
            print('\tEreverse: ', params['Ereverse'])
            print('\tomega: ', params['omega'])
            print('\tphase: ', params['phase'])
            print('\tamplitude: ', params['dE'])
            print('\tdc scan rate: ', params['v'])
            print('\ttemperature: ', params['T'])
            print('\telectrode area: ', params['a'])
            print('\tc_inf: ', params['c_inf'])
            print('\tRu: ', params['Ru'])
            print('\tCdl: ', params['Cdl'])
            print('\tE01: ', params['E01'])
            print('\tE02: ', params['E02'])
            print('\tE11: ', params['E11'])
            print('\tE12: ', params['E12'])
            print('\tE21: ', params['E21'])
            print('\tE22: ', params['E22'])
            print('\tk01: ', params['k01'])
            print('\tk02: ', params['k02'])
            print('\tk11: ', params['k11'])
            print('\tk12: ', params['k12'])
            print('\tk21: ', params['k21'])
            print('\tk22: ', params['k22'])
            print('\tGamma: ', params['Gamma'])
        except NameError as e:
            print('NameError: ', e.value)

        self.dim_params = params

        self.E0, self.T0, self.L0, self.I0 = self._calculate_characteristic_values(
        )

        self.params = dict()
        self.params['Estart'] = self.dim_params['Estart'] / self.E0
        self.params['Ereverse'] = self.dim_params['Ereverse'] / self.E0
        self.params['omega'] = 2 * pi * self.dim_params['omega'] * self.T0
        self.params['phase'] = self.dim_params['phase']
        self.params['dE'] = self.dim_params['dE'] / self.E0

        self.params['k01'] = self.dim_params['k01'] * self.T0
        self.params['k02'] = self.dim_params['k02'] * self.T0
        self.params['k11'] = self.dim_params['k11'] * self.T0
        self.params['k12'] = self.dim_params['k12'] * self.T0
        self.params['k21'] = self.dim_params['k21'] * self.T0
        self.params['k22'] = self.dim_params['k22'] * self.T0
        self.params['E01'] = self.dim_params['E01'] / self.E0
        self.params['E02'] = self.dim_params['E02'] / self.E0
        self.params['E11'] = self.dim_params['E11'] / self.E0
        self.params['E12'] = self.dim_params['E12'] / self.E0
        self.params['E21'] = self.dim_params['E21'] / self.E0
        self.params['E22'] = self.dim_params['E22'] / self.E0

        self.params['alpha1'] = self.dim_params['alpha1']
        self.params['alpha2'] = self.dim_params['alpha2']
        self.params['alpha11'] = self.dim_params['alpha11']
        self.params['alpha12'] = self.dim_params['alpha12']
        self.params['alpha21'] = self.dim_params['alpha21']
        self.params['alpha22'] = self.dim_params['alpha22']

        self.params['Ru'] = self.dim_params['Ru'] * abs(self.I0) / self.E0
        self.params['Cdl'] = self.dim_params['Cdl'] * \
            self.dim_params['a'] * self.E0 / (abs(self.I0) * self.T0)

        self.params['gamma'] = 1.0

    def dimensionalise(self, value, name):
        if name[0:1] == 'E':
            return value * self.E0
        elif name[0:1] == 'k':
            return value / self.T0
        elif name == 'omega':
            return value / (2 * pi * self.T0)
        elif name == 'phase':
            return value
        elif name == 'dE':
            return value * self.E0
        elif name == 'Ru':
            return value * self.E0 / self.I0
        elif name == 'Cdl':
            return value * self.I0 * self.T0 / (self.dim_params['a'] * self.E0)
        elif name == 'gamma':
            return value
        else:
            return NaN

    def non_dimensionalise(self, value, name):
        if name[0:1] == 'E':
            return value / self.E0
        elif name[0:1] == 'k':
            return value * self.T0
        elif name == 'omega':
            return value * (2 * pi * self.T0)
        elif name == 'phase':
            return value
        elif name == 'dE':
            return value / self.E0
        elif name == 'Ru':
            return value / self.E0 * self.I0
        elif name == 'Cdl':
            return value / self.I0 / self.T0 * (self.dim_params['a'] * self.E0)
        elif name == 'gamma':
            return value
        else:
            return NaN

    def suggest_times(self):
        final_time = self.params['Ereverse'] - self.params['Estart']
        return np.linspace(0, final_time, 1000, dtype='double')

    def simulate(self, times):
        times = np.asarray(times, dtype='double')
        current = np.empty_like(times)
        seq_electron_transfer3_explicit(self.params, current, times)
        return current, times

    def set_params_from_vector(self, vector, names):
        for value, name in zip(vector, names):
            self.params[name] = value
            self.dim_params[name] = self.dimensionalise(value, name)

    def get_params_from_vector(self, names):
        vector = np.zeros(len(names))
        for i in range(len(names)):
            vector[i] = self.params[names[i]]
        return vector

    def _calculate_characteristic_values(self):

        v = self.dim_params['v']
        T = self.dim_params['T']
        a = self.dim_params['a']

        # Faraday constant (C mol-1)
        F = 96485.3328959
        # gas constant (J K-1 mol-1)
        R = 8.314459848

        E_0 = R * T / F
        T_0 = abs(E_0 / v)

        Gamma = self.dim_params['Gamma']
        L_0 = 1.0
        I_0 = F * a * Gamma / T_0

        return E_0, T_0, L_0, I_0


class PintsModelAdaptor(pints.ForwardModel):

    def __init__(self, ec_model, names):
        self.ec_model = ec_model
        self.names = names
        self.current = pints_vector()

    def dimension(self):
        return len(self.names)

    def simulate(self, parameters, times):
        self.ec_model.set_params_from_vector(parameters, self.names)
        self.ec_model.simulate(use_times=times, use_current=self.current)
        return self.current
