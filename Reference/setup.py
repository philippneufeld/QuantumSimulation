# Philipp Neufeld, 2021-2022

from distutils.core import setup

requirements = [
    'numpy',
    'matplotlib',
    'scipy',
    'qutip',
    'tqdm',
]

setup(name='QNOSE Simulation',
      version='1.0',
      description='Reference inplementation for the QNOSE nitric oxide simulation',
      author='Philipp Neufeld',
      author_email='pneufeld@physik.pi5.uni-stuttgart.de',
      install_requires=requirements,
    )
