# HmtDB  

Open source version of the Human Mitochondrial Database [HmtDB](https://www.hmtdb.uniba.it).  

This version exactly reflects the live version available at [https://www.hmtdb.uniba.it](https://www.hmtdb.uniba.it), with a couple exceptions:  

* files containing passwords and API keys are not present due to obvious security reasons  
* the database only contains a few sequences due to available space limits.  


## Installation  

* Install the virtual environment: `virtualenv -p python3.6 venv`  
* Activate the virtual environment: `source venv/bin/activate`  
* Install required modules: `pip install -r requirements.txt`  
* Create the DB: `export FLASK_APP=app:app; flask create-db`  
* Update the DB: `flask update-db`  
* Migrate the DB: `flask migrate-db`  

## Running the DB  
`python run.py`  

When finished, deactivate the virtual environment: `deactivate`.  


