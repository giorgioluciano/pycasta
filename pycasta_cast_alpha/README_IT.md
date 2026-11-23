# PyCASTa CAST research alpha 0.1.0

Pacchetto consolidato, isolato dalla pipeline storica. Non e' una replica
certificata di CAST 1998 e NON e' ancora un comando PDB -> risultati autonomo.
Riutilizza l'export del backend CGAL esatto gia' installato. Nessuna nuova
compilazione; nessun fallback SciPy Delaunay; nessuna modifica ai sorgenti pycasta.

## Avvio Windows
Estrarre la cartella pycasta_cast_alpha in C:\github\pycasta.
Attivare l'ambiente Conda pycasta. Dalla root del progetto:

python pycasta_cast_alpha\launch.py --run-dir "C:\github\pycasta\output_isolated_tests\cast_exact_1a6w_20261117_193336_556949" --samples 128

Il launcher cerca ESATTAMENTE un cast_complex_all_dims.json e manifest.json,
e facoltativamente un topology.json. Se trova ambiguita' si ferma, non sceglie
arbitrariamente l'ultimo. Per altri percorsi usare cast_alpha.py --help.

Test:
cd pycasta_cast_alpha
python -m unittest -v test_cast_alpha

Dipendenze: numpy e scipy. Nessun download richiesto se gia' installate.

## Cosa esegue davvero
1. Controlla schema, alpha closure, celle, adiacenze opposte alle facce.
2. Ricalcola raggiungibilita' dell'infinito dai successori esportati.
3. Ricostruisce regioni usando facce condivise NON appartenenti al complesso.
4. Costruisce chiusura, bordo grezzo, interno meno alpha complex e bocche.
5. Confronta regioni, dual set, bordi e bocche con topology.json se presente,
   indipendentemente dai numeri assegnati alle regioni.
6. Calcola volume del supporto tetraedrico e centroide del supporto.
7. Con --samples >=2 stima il volume del supporto esterno a TUTTE le sfere
   del manifest, senza limitarsi alle quattro sfere del tetraedro.
8. Produce JSON, CSV, rapporto HTML offline e supporti OBJ per regioni/bocche.
9. Conserva contatti non-manifold nel JSON; non corregge o nasconde warning.

Il codice ha rivelato un dettaglio importante: collegare tutte le celle solo
per faccia condivisa NON riproduce le regioni del backend. Escludere le facce
alpha dall'adiacenza ricostruisce invece tutte le 14 regioni e le 69 bocche.
Questo e' confronto con il nostro backend, non con l'eseguibile CAST originale.

## Interpretazione delle misure
support_volume_A3 NON e' il volume fisico. free_volume_mc_A3 e' una stima
numerica del volume libero nella partizione ereditata. Non e' una misura
analitica CAST, non una validazione del modello atomico o del discrete flow.
Il ranking e' provvisorio: volume MC se attivo, altrimenti supporto geometrico.
Non si confrontano direttamente questi rank con quelli della pipeline storica.
Il centroide del supporto non e' garantito interno allo spazio libero.
L'errore MC riportato e' un errore standard di campionamento, non un intervallo
garantito: zero varianza osservata non prova esattezza. Aumentare --samples
per stabilizzare tasche piccole e rank vicini, non per cambiare i parametri CAST.
I supporti OBJ sono triangoli piani diagnostici, non superfici molecolari;
l'orientamento delle facce non e' garantito.

## Provenienza e sicurezza
Output nuovo in alpha_results/cast_alpha_TIMESTAMP. Input invariati.
JSON contiene hash SHA256 degli input, versione Python, kernel, alpha e manifest.
Si assume che manifest e complesso appartengano allo stesso run; indice e alpha
sono verificati, ma il vecchio export non contiene un hash coordinate comune:
non e' possibile provare crittograficamente la corrispondenza retroattivamente.
Raggi e preprocessing sono ereditati senza reinterpretarli o gonfiarli.
Predicati ambigui rilevati vengono rifiutati; nessuna perturbazione simbolica.
Controlli dei link dei vertici non inclusi: conservare validatore precedente.
Non esporta residui/atomi di parete con assegnazione CAST certificata.

## Test effettivamente eseguiti
9 test unitari: chiusura alpha, schema adiacenze, raggiungibilita', degenerazioni,
algebra di un tetraedro, volume 1/6, appartenenza a sfere, seed riproducibile.
Caso allegato 1a6w: ricostruite 14 regioni, 69 bocche; dual set, bordi e bocche
uguali a topology.json come insiemi. Esecuzione completa con 32 campioni/tetra.
I test sintetici controllano questi invarianti, NON coprono ancora tutti i casi
geometrici originali (tunnel, degenerazioni simboliche, misure analitiche).
Report di riferimento inclusi. I JSON originali e i dati molecolari non inclusi.

## Cosa manca per la replica originale
API verificata PDB -> backend corrente; modello atomico originale documentato;
audit indipendente delle regole esatte di flow e raggruppamento;
perturbazione simbolica; aree/volumi analitici molecolari, aperture e perimetri;
identificazione atomi delle pareti e dei bordi; confronto con riferimento CAST.
Il sorgente core.cpp allegato e' una versione precedente, non quello del backend
esatto che ha prodotto i JSON: questo pacchetto non lo sostituisce.
