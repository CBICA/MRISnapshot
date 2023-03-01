************
Installation
************

MRISnapshot can be installed using `pip`. We suggest that users install the package in a Conda environment or Python3 virtual environment:

**1. Installation in a Conda environment using pip**

    .. code-block:: console

        conda create -n mrisnapshot python=3.8
        conda activate mrisnapshot
        conda install pip
        
        pip install mrisnapshot

**2. Installation in a Python3 virtual environment using pip**

    .. code-block:: console
        
        virtualenv -p python3 mrisnapshot
        source mrisnapshot/bin/activate
        
        pip install mrisnapshot



**3. Installation in a Conda environment from Github repository**

    .. code-block:: console

        git clone https://github.com/CBICA/MRISnapshot
        cd MRISnapshot

        conda create -n mrisnapshot python=3.8
        conda activate mrisnapshot
        conda install pip
        
        pip install -r requirements.txt
        pip install .


**4. Check the installation**

    .. code-block:: console

        mrisnapshot_prep_data -h
        mrisnapshot_create_report -h


