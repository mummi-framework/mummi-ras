'''
Mock methods to run AA analysis in non-interactive mode.
'''

class MockProcess():

    def __init__(self):
        self.counter = 0
        self.response = None  # initialize the process to None

    def poll(self):
        self.response = None

        self.counter += 1

        return self.response


class MockAAsimulation():

    def __init__(self):
        self.curIdx = 0

    def get_curIdx(self, trajname):
        self.trajname = trajname
        self.curIdx = int(self.trajname.split(".")[1])

        return self.curIdx

    def mdcrd_checkstate(self):
        pass

    def get_islocal(self):
        return False

    def stop(self):
        pass
