************
Installation
************

MRISnapshot can be installed using `pip`. We suggest that users install the package in a Python3 virtual environment (or alternatively in a Conda environment). Installation steps are:




1. Clone the project.

.. code-block:: console

    git clone https://github.com/CBICA/MRISnapshot

2. Go to the root directory of the project.

.. code-block:: console

    cd MRISnapshot


3. Create and activate a Python 3 virtual environment

.. code-block:: console

    virtualenv -p python3 mrisnapshot
    source mrisnapshot/bin/activate
    
.. note::
    Alternatively, users can also create a conda environment

    .. code-block:: console

        conda create -n myenv python=3.8
        conda activate mrisnapshot
        conda install pip
    
4. Install the package

We recommend installation of the latest version of the package directly from pypi:

.. code-block:: console
    
    pip install mrisnapshot
    

Alternatively, the package can be install from the downloaded repository:

.. code-block:: console

    pip install -r requirements.txt
    pip install .


5. Check the installation

.. code-block:: console

    mrisnapshot_prep_data -h
    
    mrisnapshot_create_report -h


