# Updating protocol  

## Dependencies installation  

Seguire quanto descritto nel file INSTALL.txt.  

## Updating HmtDB  

### Ricerca nuovi genomi  

**Script `genome_search.py`**  

Lo script si occupa di ricercare in GenBank i nuovi genomi mitocondriali umani pubblicati a partire dalla data dell'ultimo aggiornamento del db. È comunque possibile specificare una data da utilizzare per la ricerca dei nuovi genomi, usando l'opzione `-f/--from_date 10_2017`.  
Un argomento richiesto è la fase della ricerca: 
* inserendo `1` verranno cercati i nuovi genomi, parsate le loro informazioni che saranno salvate in due tabelle csv (GenBank e PubMed), nonché salvati i singoli file GenBank e Fasta per ogni nuovo genoma trovato;  
* a questo punto va modificato il csv GenBank per inserire nell'ultima colonna `Y` o `N` se il genoma è sano o diseased, quindi rilanciare lo script con l'argomento `2`, in modo che venga prodotto anche il multifasta in cui saranno già inseriti i nuovi HmtDB IDs dei genomi trovati.  

I file `log_genome_search_[ph1/ph2]_mm_yyyy.log` contengono tutti i dettagli dell'update eseguito.  
Usage:  
```
python genome_search.py 1  
python genome_search.py 2  
```  
**N.B.: La fase 2 dello script va lanciata solo dopo aver concluso la parte di ML eseguita dallo script `ml_predict.py`.**  

#### Genome classification with ML  

**Script `ml_predict.py`**  

Questo script si occupa di distinguere genomi healthy e patient usando un approccio di machine learning. L'informazione viene recuperata dall'abstract della pubblicazione in cui è stato riportato ogni genoma; l'algoritmo viene istruito sulla base dei genomi già presenti nel db, e impara così una serie di parole chiave che distinguono gli abstract dei genomi healthy da quelli dei patient, sulla base delle quali riesce poi a classificare i nuovi genomi.  
Se lanciato senza opzioni, lo script riproduce automaticamente il file di training sulla base delle informazioni ottenute dall'ultimo update.  
Va usato dopo la fase 1 dello script `genome_search.py` in modo da avere l'informazione necessaria per poter poi lanciare la fase 2.  

Il file `log_ml_predict_mm_yyyy.log` contiene tutti i dettagli.  
Usage:  
```
python ml_predict.py
``` 
Viene prodotto il file `model_output.csv`, che contiene i PubmedId analizzati, la predizione (0 per genomi normali, 1 per genomi pazienti) e la probabilità di questa predizione. Questi risultati andrebbero controllati, almeno quelli in cui quest'ultimo valore risulta minore di ~0.7. Fatto questo, bisogna riportare i risultati nel csv GenBank prodotto nel primo step, inserendo Y dove i genomi sono sani e N altrimenti.  


### Allineamento automatizzato dei nuovi genomi  

**Script `align_seqs.py`**  

Lo script effettua l'allineamento automatizzato tra il multifasta con i nuovi genomi creato nello step precedente ed i genomi già allineati presenti nel db. L'allineamento è eseguito da Mafft, sulla base di un template rappresentato dai genomi di HmtDB già allineati, in modo da preservare l'allineamento locale di alcune zone particolari (stretch ecc.).  
L'allineamento creato viene poi caricato sul db.  
Usage:  
```
python align_seqs.py
```  

### Caricamento dati sul db  

**Script `load_entries.py`**  

Questo script carica nel db le nuove tabelle prodotte dalle fasi 1 e 2 dello script `genome_search.py`, in modo da poter poi lanciare il calcolo delle variability.  
Usage:  
```
python load_entries.py
```

### Nucleotide variability  

**Script `calc_sitevar.py`**  

Lo script crea di default un file csv chiamato `sitevar_healthy` o `sitevar_patient`; se lanciato senza opzioni calcola la variabilità totale dei pazienti, con `-H` (`--healthy`) quella totale dei sani.  
Per effettuare il calcolo solo su uno specifico continente, aggiungere il flag `-C` (`--continent`) seguito da AF, AM, AS, EU, OC.  

**Script `nt_var.py`**  

Si tratta di un wrapper dello script precedente, che calcola le variabilità di sani e pazienti, totale e continente-specifica, per poi caricare i dati sul db.  
Va lanciato su Recas, in quanto ci mette un po'.  

Per il classico aggiornamento generale di HmtDB, va lanciato come `python nt_var.py -c`; per il calcolo delle variability dei completi (per l'aggiornamento di HmtVar) va usato come `python nt_var.py -c -l -M`.  

### Allele frequencies  

**Script `allele_freqs.py`**  

Questo script si occupa di calcolare le frequenze alleliche (dato che a quanto pare la porzione dell'`nt_var.py` che anche se dovrebbe farlo non funziona a dovere), le salva in specifici file csv nella cartella `all_freqs` e le carica nel db.  
Per lanciarlo sui sani:  
```
python allele_freqs.py -H
```
Per lanciarlo sui pazienti basta eliminare il flag `-H`.  
Per calcolare le frequenze alleliche sui genomi completi ed evitare di caricare questi dati sul db, basta lanciare lo script usando:  
```
python allele_freqs.py -C -l
```
ed aggiungere il flag -H quando necessario.  

### Aminoacid variability  

**Script `aa_var.py`**  

Calcola la variabilità aminoacidica e carica i dati sul db.  
Va lanciato in due fasi successive: durante la prima viene prodotto il file `aa_var_temp.sh`, che deve poi essere eseguito usando qsub. Nella fase 2, i risultati di quest'ultimo processo vengono caricati sul db. Questo iter va ripetuto per i sani e per i pazienti.  
In sostanza, il protocollo sarebbe il seguente:  

``` 
python aa_var.py -s 1     # genera aa_var_temp.sh per i pazienti (va lanciato interattivamente, NON CON QSUB, tanto va veloce
qsub aa_var_temp.sh       # calcola aa var per i pazienti 
python aa_var.py -s 2     # carica i dati sul db 
python aa_var.py -s 1 -H  # genera aa_var_temp.sh per i sani (va lanciato interattivamente, NON CON QSUB, tanto va veloce
qsub aa_var_temp.sh       # calcola aa var per i sani 
python aa_var.py -s 2 -H  # carica i dati sul db 
```

### Zip data  

**Script `zip_data.py`**  

Salva i multiallineamenti ed i dati di variabilità in file zippati scaricabili da HmtDB.  

```
python zip_data.py
```

### Calculate statistics  

**Script `calc_stats.py`**  

Calcola le statistiche visualizzate poi sul sito di HmtDB.  

```
python calc_stats.py
```

