Installing **dnora**
=============================================

1. Using pip (recommended)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The easiest way to install dnora is to use pip. It is recommended to create a new conda environment to work in:

.. code-block:: bash

   $ conda create -n dnora python=3.10
   $ conda activate dnora
   $ python -m pip install dnora

2. Cloning the GitHub repository (for development)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
If you want to work on the code (or use some other branch than main), you can pull take the original code:

1. Pull the code from GihHub

.. code-block:: bash

   $ git clone https://github.com/MET-OM/dnora.git
   $ cd dnora/

2. Create an environment with the required dependencies and install dnora as an editable package

.. code-block:: bash

  $ conda config --add channels conda-forge
  $ conda env create -f environment.yml
  $ conda activate dnora2
  $ pip install -e .
  
3. Installing with fimex
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Installing dnora with fimex might be a bit tricky. Fimex is not in pip, and dnora is not in conda. You need to install fimex first, and dnora second, otherwise you run into dependecy issues:

.. code-block:: bash

   $ conda create -n dnora python=3.10 fimex=1.8.1
   $ conda activate dnora
   $ python -m pip install dnora

Depending on the platform you might need to do something like:

.. code-block:: bash

   $ mamba create -n dnora python=3.12 fimex -c conda-forge
   $ mamba activate dnora
   $ python -m pip install dnora
