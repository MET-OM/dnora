Installing **DNORA**
=============================================

Quick Installation: 
The easiest way to install dnora is to use pip:

.. code-block:: bash

   $ pip install dnora

Alternatively, you can clone the GitHub repository:

1. Install anaconda3 or miniconda3
2. Clone dnora:

.. code-block:: bash

   $ git clone https://github.com/MET-OM/dnora.git
   $ cd dnora/

3. Create an environment with the required dependencies and install dnora

.. code-block:: bash

  $ conda config --add channels conda-forge
  $ conda env create -f environment.yml
  $ conda activate dnora2
  $ pip install .
  
To update the environment using a new environment.yml, run:

.. code-block:: bash

   $ conda env update --file environment.yml --prune

For JupyterHub platforms such as WEKEO, please use:

.. code-block:: bash

   $ conda create -n dnora python=3.10 fimex=1.8.1
   $ conda activate dnora
   $ python -m pip install dnora