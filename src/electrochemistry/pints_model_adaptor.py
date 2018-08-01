class PintsModelAdaptor(pints.ForwardModel):

    def __init__(self, ec_model, names):
        self.ec_model = ec_model
        self.names = names

    def n_parameters(self):
        return len(self.names)

    def n_outputs(self):
        return 1

    def simulate(self, parameters, times):
        self.ec_model.set_params_from_vector(parameters, self.names)
        return self.ec_model.simulate(times)


