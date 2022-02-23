# Philipp Neufeld, 2021-2022
# pip install ARC-Alkali-Rydberg-Calculator .

import matplotlib.pyplot as plt
import numpy as np
import arc

if __name__ == "__main__":
    
    atom = arc.Caesium()
    calc = arc.StarkMap(atom)
    calc.defineBasis(28, 0, 0.5, 0.5, 23, 32, 20)
    calc.diagonalise(np.linspace(00.,60000,600))
    calc.plotLevelDiagram()
    calc.showPlot()

