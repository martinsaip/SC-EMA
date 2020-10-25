__author__ = 'aydin'

myStack = {}


def put(value = None):
    if value == None:
        print("empty value, nothing to do")
        return
    else:
        myStack['key_'+str(len(myStack))] = value

def showStack():
    print(myStack)

def executeStack(key = None):
    eval(myStack[key])

class crystalStack:
    def __init__(self):
        self.stack = {}

    def put(self, value = None):
        if value == None:
            print("empty value, nothing to do")
            return
        else:
            #TODO: something to do with self.stack
            print(eval(value))

    def get(self):
        return self.stack
