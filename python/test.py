import dnsbox as db
import numpy as np
db.init()
db.setTimeStep(0.001)
db.saveState()

for n in range(10000):
    db.timeStep(10)
    Ekin = db.Ekinetic()
    print ('t = ', db.t)
    print ('Ekin =', Ekin)

    if n == 5:
        db.saveState()

# db.checkError()