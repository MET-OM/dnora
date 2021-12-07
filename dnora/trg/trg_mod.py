from copy import copy

class TrGrid:
    def __init__(self, name='AnonymousTrianGrid'):
        self._name = copy(name)
        return


    def name(self):
        return self._name
