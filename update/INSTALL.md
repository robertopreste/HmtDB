# Istruzioni per l'installazione del protocollo di update  

N.B.: per motivi storici l'update di HmtDB richiede Python 2.7. La maniera migliore per 
soddisfare questa dipendenza è installarla usando conda con `conda create --name py2 python=2.7 virtualenv`. 
Creare una virtualenv con Python 2.7: `virtualenv -p python2.7 venv_upd`. 

Tutti i seguenti step vanno eseguiti ovviamente dopo aver attivato la virtualenv usando `source venv_upd/bin/activate` ed aver cambiato working directory in db_update `cd db_update`.  

## MAFFT (multiallineamento)  
Info disponibili a [questo link](https://mafft.cbrc.jp/alignment/software/installation_without_root.html)  
* `wget https://mafft.cbrc.jp/alignment/software/mafft-7.394-without-extensions-src.tgz`  
* `tar -zxvf mafft-7.394-without-extensions-src.tgz`  
* `mv mafft-7.394-without-extensions/ mafft-linux`  
* `cd mafft-linux/core`  
* `nano Makefile`  
* modificare: `PREFIX = /this/local/folder/mafft-linux`  
* `make clean`  
* `make && make install`  
* `cd ../../`

## Bioinf_pkg e Pyrex (variabilità nucleotidica)  
* decomprimere il file Bioinf_pkg.tar.gz: `tar zxvf Bioinf_pkg.tar.gz`  
* OPPURE scaricare la repo da [https://github.com/robertopreste/Bioinf_pkg](https://github.com/robertopreste/Bioinf_pkg)  
* OPPURE clonare direttamente la repo: `git clone git@github.com:robertopreste/Bioinf_pkg.git`  
* `cd Bioinf_pkg/`  
* Scaricare il pacchetto Pyrex: `wget http://www.cosc.canterbury.ac.nz/greg.ewing/python/Pyrex/Pyrex-0.9.9.tar.gz`  
* `tar -zxvf Pyrex-0.9.9.tar.gz`  
* `cd Pyrex-0.9.9/`  
* `python setup.py install`  
* `cd ../`  
* `python setup.py build`  
* `python setup.py install`  
* `cd ../`  
* `pip install -r ../requirements.txt`  
* `conda install pandas`  

## Mitvar (variabilità aminoacidica)  

* `cd mitvar`  
* rimuovere l'eseguibile se già presente: `rm mitvarprot`  
* compilare l'eseguibile: `g++ mitvarprot.cpp -o mitvarprot`  
* `cd ../`

## EMBOSS (transeq e revseq necessari a mitvarprot)  

* `wget ftp://emboss.open-bio.org/pub/EMBOSS/emboss-latest.tar.gz`  
* `tar zxvf emboss-latest.tar.gz`  
* `mv EMBOSS-6.6.0/ emboss`  
* `cd emboss`  
* `./configure --prefix=<this/local/folder/emboss> --without-x`  
* `make && make install`  
* `cd ../`  

