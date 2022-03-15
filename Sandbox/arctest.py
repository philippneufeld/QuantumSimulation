# Philipp Neufeld, 2021-2022
# pip install ARC-Alkali-Rydberg-Calculator .

import matplotlib.pyplot as plt
import numpy as np
import arc

if __name__ == "__main__":
    
    atom = arc.Hydrogen()
    
    # levels = arc.LevelPlot(atom)
    # levels.makeLevels(23,32,0,20)
    # levels.drawLevels()
    # levels.showPlot()

    calc = arc.StarkMap(atom)
    calc.defineBasis(28, 0, 0.5, 0.5, 23, 32, 20)

    calc.eFieldCouplingSaved = arc.calculations_atom_single._EFieldCoupling()
    print(atom.getRadialMatrixElement(28, 2, 2.5, 28, 3, 2.5) * 0.5291772083e-10)
    print(calc._eFieldCouplingDivE(28, 2, 2.5, 0.5, 28, 3, 2.5, 0.5))

    # n, l, j = 28, 2, 2.5
    # a1,b1 = atom.radialWavefunction(l,0.5,j, atom.getEnergy(n, l, j)/27.211, atom.alphaC**(1/3.0), 2.0*n*(n+15.0), 0.001)
    # n, l, j = 28, 3, 2.5
    # a2,b2 = atom.radialWavefunction(l,0.5,j, atom.getEnergy(n, l, j)/27.211, atom.alphaC**(1/3.0), 2.0*n*(n+15.0), 0.001)
    # plt.plot(a1, a1**2*b1**2)
    # plt.plot(a2, a2**2*b2**2)
    # plt.show()

    print(calc.basisStates)
    for row in calc.mat2:
        print("\t".join(f"{el:.3e}" for el in row * 1e9 * 6.62607004e-34))

    calc.diagonalise(np.linspace(00.,60000,600))
    calc.plotLevelDiagram()
    calc.showPlot()

