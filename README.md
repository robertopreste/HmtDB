# HmtDB  

The Human Mitochondrial Database [HmtDB](https://www.hmtdb.uniba.it).  


## Installation  

Only the first time the database is set up:  

1. Create a new virtual environment (Python 3): `virtualenv -p python3.6 venv`  
2. Activate the virtual environment: `source venv/bin/activate`  
3. Install required modules: `pip install -r requirements.txt`  
4. Create the DB: `export FLASK_APP=app:app; flask create-db`  
5. Update the DB: `flask update-db`  
6. Migrate the DB: `flask migrate-db`  

When finished, deactivate the virtual environment: `deactivate`.  

## Updates  

Please refer to the [/update/INSTALL.md](/update/INSTALL.md) and [/update/README.md](/update/README.md) files for details about the updating protocol for HmtDB.  
After the updating procedure is complete, please run `flask migrate-db`.  

## HmtDB instance  

HmtDB is served using [gunicorn](https://gunicorn.org); please ask your system admin for help about this.  
