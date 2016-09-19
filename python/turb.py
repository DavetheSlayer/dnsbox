import dnsbox as db
import numpy as np
db.init()

for n in range(100):
    db.stats()            # Compute and save flow statistics
    db.setTimeStep()      # Set time step
    if n % 100 == 0:
        db.saveState()    # Save state every 1000 steps
    db.timeStep(10)

db.setTimeStep()      # Save final step
# db.checkError()