## Algoritmi di programmazione matematica per il Max Edge Weighted Clique problem with Multiple Choice constraints.

#### Questo progetto rappresenta il lavoro svolto per la tesi di laurea magistrale in Scenze e tecnologie dell'informazione.

Molti problemi rilevanti in ambito scientifico e tecnologico possono essere descritti e risolti con successo utilizzando modelli basati su grafi. Ad esempio, in biologia computazionale sorge frequentemente il problema conosciuto come side-chain placement problem che consiste nel trovare la conformazione lineare di una proteina, in cui l'energia di legame tra gli amminoacidi sia minima e che, al contempo, assicuri un'adeguata stabilità strutturale. In modo simile, nell'ambito della management science, per garantire l'interoperabilità dell'insieme di reparti di una grande impresa è necessario scegliere il protocollo di interazione che massimizzi l'efficienza complessiva. Anche nell'ambito del progetto di reti di telecomunicazione sorge spesso l'esigenza di disporre i nodi che compongono la dorsale in modo da minimizzare i tempi di latenza o massimizzare la banda passante. Tutti questi problemi hanno una struttura comune, riconducibile alla ricerca di una clique di peso massimo in un grafo opportuno. In questa tesi viene a#rontato il Max Edge Weighted Clique problem with Multiple choice Contraints (MEWCMC) che, dati un grafo pesato sia sui vertici che sui lati, e una partizione dell'insieme dei vertici, consiste nell'individuare una clique di peso massimo selezionando un vertice per ogni classe della partizione. Il MEWCMC è un problema di ottimizzazione su grafo molto generico e, come molti altri problemi simili, è N P − difficile. Attualmente non sono noti algoritmi in grado di risolvere in tempo polinomiale questa classe di problemi: per ottenere metodi di risoluzione efficaci, è necessario studiare il problema in modo approfondito. Nella tesi vengono proposti modelli matematici e algoritmi per la risoluzione esatta del problema, basati su tecniche di enumerazione implicita. In particolare è stato dapprima costruito un modello di programmazione lineare intera, un modello con funzione obiettivo quadratica e uno di programmazione semidefinita ottenuto tramite riformulazione del modello quadratico. Un'analisi sperimentale ha evidenziato che i solutori commerciali, utilizzando questi modelli, sono in grado di risolvere solo istanze di dimensioni ridotte. Quindi, per superare le limitazioni dei solutori generici è stato sviluppato e implementato un algoritmo euristico di Tabu Search, in grado di determinare delle soluzioni di buona qualità, in tempi molto rapidi. Infine è stato studiato e implementato un algoritmo ad-hoc, che sfrutta euristiche e tecniche di rilassamento per esplorare efficientemente l'albero di decisione. Il primo ingrediente dell'algoritmo è un rilassamento che utilizza tecniche di programmazione semidefinita per ottenere una stima per eccesso del valore della soluzione ottima. Il secondo ingrediente è un bound che sfrutta unicamente le proprietà combinatorie del problema e può essere calcolato in modo molto efficiente. Tale bound è utilizzato anche per ridurre le dimensioni iniziali del problema,eliminando a priori soluzioni sub-ottime. Il terzo ingrediente è un insieme di politiche di branching volte a migliorare l'efficienza nell'esplorazione dell'albero di ricerca. L'ultimo ingrediente è un'euristica di rounding, che utilizza le soluzioni fornite dal rilassamento per individuare altre soluzioni ammissibili. Poiché la tecnica di rilassamento gioca un ruolo cruciale nell'efficienza dell'algoritmo ad-hoc, è stata dedicata particolare attenzione alla possibilità di introdurre vincoli aggiuntivi, in grado di migliorare il bound senza aumentare significativamente il carico computazionale. L'implementazione ha richiesto lo sviluppo di codice specifico con l'integrazione di solutori esterni. Infatti, la valutazione del rilassamento basato su programmazione semidefinita richiede l'impiego di algoritmi specifici, fra cui il metodo del cammino centrale, la cui efficienza varia sensibilmente in funzione della scelta di alcuni parametri. 
Le prestazioni degli algoritmi implementati sono state valutate tramite un'estesa campagna sperimentale. Questa ha evidenziato che l'impiego delle tecniche di programmazione semidefinita è particolarmente efficace: permette di ottenere bound duali stretti, riducendo così il tempo necessario all'enumerazione implicita. 
I risultati ottenuti dall'algoritmo ad-hoc sono stati confrontati con quelli dal solutore commerciale ILOG Cplex 11.2 a cui è stato fornito in ingresso il modello lineare intero. L'algoritmo ad-hoc risolve istanze di dimensioni molto maggiori rispetto a Cplex, comportandosi meglio sia in termini di tempo di calcolo che in termini di qualità dei bound. È stato inoltre possibile verificare che l'euristica di Tabu Search implementata, nella maggior parte delle istanze testate, riesce ad individuare la soluzione ottima in pochi secondi. L'efficacia del rilassamento semidefinito suggerisce che altre tecniche di rilassamento, basate su programmazione conica, possano essere efficaci nella soluzione di MEWCMC così come di altri problemi che presentano una struttura simile.

-----

Aiutandomi con la presentazione disponibile nella cartella doc vi riassumo brevemente il lavoro.

Definizione del problema:
Dato un grafo partizionato, completo, pesato sia sui vertici che sui lati, il problema consiste nel determinare una clique di peso massimo. I vincoli di multiple choice impongono che venga selezionato esattamente un vertice per ogni classe.
La figura mostra un esempio molto semplice di come è calcolato il valore della clique formata dai vertici 2,5,8
Si nota che tale problema contiene come caso particolare il problema della clique di peso massimo.

Questo problema sorge in numerose applicazioni reali fra cui il Data Mining, la biologia computazionale, le reti di telecomunicazioni e il managemnet science. In particolare nella tesi abbiamo presentato 4 scenari reali in cui il problema potrebbe trovare applicazione, ma ce ne potrebbero essere molte altre.

Letteratura:
In letteratura è presente una lavoro del 2005 in cui viene presentato un branch and bound con rilassamento lagrangeano.
Il problema è molto studiato nella versione senza i vincoli di MC, sino dal 1977
Da allora sono stati presentati numerosi modelli, e algoritmi per risolverlo, sia euristiche basate su Tabu Search e Grasp (Greedy randomized adaptive search procedure) che esatti basati. Fra gli algoritmi esatti sono stati proposti dei branch and bound con piani di taglio.

In questa tesi abbiamo proposto dei nuovi modelli per il problema, degli algoritmi euristici, e un algoritmo esatto che include di bounding innovative.
Ho verificato le prestazioni degli algoritmi con una estesa campagna sperimentale e ho confrontato i risultati con la letteratura e con il solutore commerciale cplex.

Un modello matematico quadratico:
Una formulazione molto semplice pe il problema si ottiene introducendo n variavili decisionali binarie xi, una per ogni vertice.
Ogni termine della funzione obiettivo è diverso da zero se entrambe le variabili xi e xj sono a 1 ovvero entrambi i vertici associati fanno parte della clique.
I vincoli 1 di multiple choice impongono che venga selezionato esattamente un vertice per ogni classe della partizioni.

Questo modello è difficile da risolvere poiché non è convesso: i vincoli di integralità inoltre se rilassiamo tali vincoli la funzione obiettivo è convessa solo se la matrice W è semidefinita positiva.

Possimao ottenere un modello lineare intero, più facile da risolvere introducendo n^2 nuove variabili decisionali binarie y_ij che prendono il valore del prodotto xi*xj. Per imporre questa condizione è necessario introdurre 3 nuovi gruppi di vincoli.
Il rilassamento continuo del modello ottenuto formissce un upper bound per il problema.

Modello matematico con le matrici:
organizzando le variabili y in una matrice Y possiamo ottenere un altro modello di programmazione matematica basata su matrici
A questo punto riscrivo la funzione obiettivo come il prodotto di frobenius tra la matrice peso e la matrice incognita Y, rappresento i vincoli di multiple choice mediante opportune matrici.

I vincoli 8-10 impongono rispettivamente che Y abbia rango 1 e che sia semidefinita positiva. Tecnicamente questi vincoli mi sevono per imporre che Y possa essere scritta come Y=x*x'  ovvero le righe e le colonne sono tutte linearmente dipendenti e gli elementi sono binari.

Euristiche:
Trovare una soluzione ammissibile per il problema è semplice, basta selezionare un vertice per ogni classe ma trovare la clique di peso massimo è molto difficile.
Abbiamo proposto alcuni algoritmi euristici: il primo si basa esclusivamente sulle proprietà combinatorie del problema.
Per ogni vertice definiamo un contributo potenziale che tiene conto del peso del vertice, dei pesi dei lati incidenti e dei vertici appartenenti alle altre classi. Schegliendo per ogni classe il vertice con il maggior contributo potenziale otteniamo una soluzione ammissibile.

Abbiamo proposto anche un'altra euristica molto efficiente basata su Tabu Search che parte da una soluzione iniziale ed esplora di volta in volta l'intorno della soluzione costituito da tutti i possibili scambi di un vertice nella solzuione un uno della stessa classe che non appartiene. La scelta della soluzione dell'intorno viene fatta con politica best improve.
Il Tabu Search fornisce una buona soluzione di partenza per l'algoritmo esatto

Dato che i modelli presentati forniscono delle soluzioni frazionarie abbiamo implementato delle tecniche di arrotondamento per ottenere delle soluzioni ammissibili.

Branch & Bound algorithm:
L'algoritmo di sviluppato fa parte della classe degli algoritmi di branch and bound e consiste nell'esplorazione di un albero di decisione in cui i nodi rappresentano i spottoproblemi. Le tecniche di bounding formiscono delle stime per eccesso e per difetto del valore della solzione ottima di un sottoproblema e permettono di stabilire a priori quali sottoproblemi tralasciare.

Di seguito presentiamo le componenti dell'algoritmo, in particolare, le tecniche di bounding, sia duale che primale, la tecnica di branching e altre tecniche che permettono di aumentare l'efficienza.

Il rilassamento semidefinito:
Consideriamo il modello precedente basato sulle matrici. Rilassando il vincolo di rango otteniamo un problema di ottimizzazione convessa nel continuo per cui esistono numerose tecniche basate sul metodo di Newton in grado di determinare l'ottimo globale.
La soluzione rappresenta un bound duale per il problema.
In particolare il problema così formulato è un problema di programmazione semidefinita per cui esistono delle tecniche molto efficienti ovvero degli algoritmi di tipo primale-duale basati sul metodo del punto interno. 
Per risolvere il modello abbiamo utilizzato il solutore sperimentale DSDP che implementa l'algoritmo del cammino centrale per determinare la matrice soiluzione Y.

Alla soluzione frazionaria ottenuta possiamo applicare la dell'arrotondamento per ottenere una soluzione primale valida.

Bound duale combinatorio:

Utilizziamo l'idea già vista per l'euristica primale combinatoria per definire un bound duale combinatorio.

Per ogni classe consideriamo il vertice che ha il massimo contributo potenziale in modo indipendente dalle altre classi.
L'unione dei contributi delle classi fornisce una stima per eccesso del valore della soluzione ottima del problema.

Nella tesi abbiamo mostrato che questa stima è un bound duale valido mediante la tecnica della maggiorazione.
Sfruttando le proprietà combinatorie abbiamo inoltre implementato un algorirmo di preprocessing in grado di ridurre le dimensioni del problema di partenza.

Strategia di branching:
La politica di branching permette di definire come vengono generati i sottoproblemi in particolare in questa slide riportiamo un esempio della soluzione ottenuta dal rilassamento del modello semidefinito. Sulla diagonale della matrice si trovano i valori frazionari rilativi ai vertici.

Ho implementato una strategia che permette di generare due sottoproblemi figli partizionando le variabili frazionarie di una classe in due sottoinsiemi: nel primo figlio il secondo sottoinsieme viene fissato a 0 analogamente per il secondo figlio, le variabili del primo sottoinsieme vengono fissate a 0. Possiamo quindi eliminare dalla matrice le righe e le colonne relative alle variabili fissate.

Come risultato il problema viene ridotto di dimensione e quindi è più facile e conseguentemente il bound duale si stringe.

Enumerazione esplicita:
Il branching riduce progressivamente la dimensione del problema... per ogni sottoproblema  il numero di soluzione ammissibili è pari al prodotto delle variabili libere in ogni classe.
Nonostante le tecniche di bounding siano polinomiali, sono computazionalmente onerore e quindi esiste una dimensione del problema per cui risulta conveniente enumerare esplicitamente le soluzioni.

Esperimenti computazionali:
Tutti gli algoritmi sono stati scritti in C in ambiente cplex.
Per testare le prestazioni degli algoritmi è stata condotta una vasta campagna di simulazioni
abbiamo inoltre confrontato i risultati con due competitor: il solutore commerciale cplex nella sua ultima versione a cui è stato fornito il modello di PLI con vincoli di taglio e l'algoritmo lagrangeano proposto in letteratura.

Abbiamo considerato due dataset di istanze: il primo creato da noi e il secondo preso dalla letteratura.
Riportiamo solo i risultati più significativi in una forma aggregata.

Primo grafico:
In questa slide presentiamo i risultati realtivi al dataset 1  che è composto da istanze fino a 300 vertici diviso in 3 classi con diverse caratteristiche nella generazione dei pesi dei vertici/lati.

Il nostro algoritmo si è dimostrato efficiente, in un ora è stato in grado di chiudere il 77% delle istanze mentre cplex solo il 44.

Il primo grafico mostra i tempi di calcolo per ogni classe quando entrambi gli algoritmi trovano la soluzione ottima.
Si nota un ordine di grandezza di differenza in meno rispetto a cplex, puntualmente anche 2-3

Il secondo grafico considera i gap di ottimalità. E' ragionevole pensare che dando un po' di tempo al nostro algoritmo è possibile chiudere tutte le istanze al contrario. Dai risultati della classe C concludiamo che non vale la stessa cosa per cplex.

Secondo grafico:
Mostriamo ora i risultati sul dataset 2 proposto in letteratura composto da 168 istanze divise in 4 classi con diverse caratteristiche.
Il grafico mostra una nuova colonnina rosa relativa al branch and bound lagrangiano. I risultati sono in linea con il grafico precedente ovvero il nostro algoritmo è più veloce di 1-2 ordini di grandezza rispetto ai competitor.

Conclusioni:
In questo lavoro di tesi abbiamo studiato le proiprietà copmbinatorie del problema, abbiamo proposto dei modelli e studiato delle tecniche di bounding innovative.
Le prestazioni mostrano tempi di calcolo fino a 2-3 ordini di grandezza rispetto ai competitor.
Abbiamo inoltre verificato che l'euristica di TS è stata in grado di determinare la soluzione ottima nel 84% delle istanze.

Dato che il problema studiato ha come caso particolare la clique di peso massimo ci sono numerosi problemi che vi possono essere ricondotti tramite l'introduzione di vincoli aggiuntivi.

Il successo delle tecniche di bounding suggerisce l'utilizzo di tecniche vbasate sul second order conic programming, ancora fortemente sperimentali.

Qualora foste interessati al dataset delle istanze utilizzate per gli esperimenti computazionali contattatemi.
