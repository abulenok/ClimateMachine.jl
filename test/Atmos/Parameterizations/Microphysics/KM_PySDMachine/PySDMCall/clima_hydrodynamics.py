
"""
Class for utilizing ClimateMachine's hydrodynamics
"""


class ClimateMachine:

    def __init__(self, clima_fields):
        self.fields = {}
        self.fields['qv'] = clima_fields['qv']
        self.fields['thd'] = clima_fields['thd']
        self.particulator = None

    def register(self, builder):
        self.particulator = builder.particulator

    def __call__(self):
        self.particulator.env.get_predicted('qv').download(self.particulator.env.get_qv(), reshape=True)
        self.particulator.env.get_predicted('thd').download(self.particulator.env.get_thd(), reshape=True)

    def set_qv(self, qv):
        self.fields['qv'] = qv

    def set_thd(self, thd):
        self.fields['thd'] = thd
