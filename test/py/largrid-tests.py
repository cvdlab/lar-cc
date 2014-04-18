import os
def createDir(dirpath):
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)

createDir('test/py/largrid/')
