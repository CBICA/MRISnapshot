************
Installation
************

MRISnapshot can be installed using `pip`. We suggest that users install the package in a Python3 virtual environment. Installation steps are:

1. Clone the project.

.. code-block:: console

    git clone https://github.com/CBICA/MRISnapshot

2. Go to the root directory of the project.

.. code-block:: console

    cd MRISnapshot


3. Create and activate a Python 3 virtual environment

.. code-block:: console

    virtualenv -p python3 MRISnapshot
    source MRISnapshot/bin/activate

4. Install package requirements

.. code-block:: console

    pip install -r requirements.txt

5. Install the package

.. code-block:: console

    pip install .

5. Check the installation

.. code-block:: console

    mrisnapshot_prep_data -h
    
    mrisnapshot_create_report -h


