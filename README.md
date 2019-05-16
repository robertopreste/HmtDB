# HmtDB  

Open source version of the Human Mitochondrial Database [HmtDB](https://www.hmtdb.uniba.it).  

This version exactly reflects the live version available at [https://www.hmtdb.uniba.it](https://www.hmtdb.uniba.it), with a couple exceptions:  

* files containing passwords and API keys are not present due to obvious security reasons  
* the database only contains a few sequences due to available space limits.  


## Installation  

* Install the virtual environment: `virtualenv -p python2.7 venv`  
* Activate the virtual environment: `source venv/bin/activate`  
* Install required modules: `pip install -r requirements.txt`  

## Creating the database  

```bash
mysql -u root -p  # root password
```

Create a new user for HmtDB:  

```mysql
USE mysql; 
CREATE USER 'hmtdb_admin'@'localhost' IDENTIFIED BY 'password';
GRANT ALL PRIVILEGES ON *.* TO 'hmtdb_admin'@'localhost';
FLUSH PRIVILEGES; 
```

Exit MySQL (using `\q`) and enter back using the new credentials:  

```bash
mysql -u hmtdb_admin -p password
```

Create the database:  

```mysql
CREATE DATABASE HmtDB; 
```

## Migration and upgrade  

Export the Flask app using `export FLASK_APP=app:app`, then:  

* create the database: `flask create-db`  
* or open a Python interpreter and issue:  
```python
from app import db 
db.drop_all()
db.create_all()
```
* load all the tables into the db: `flask update-db`  

TODO: update tables and load them on the DB. 

## Running the DB  
`python run.py`  

When finished, deactivate the virtual environment: `deactivate`.  


