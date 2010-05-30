#!/bin/bash

#parametri:
# numero iterazioni tabu
# file output


#Simulazione su tutte le istanze

path_istanze=./istanze
fileout="./output_simulazione_dsdp"
fileout=$2 
echo "Simulazione dsdp" > $fileout;



for classe in {a,b,c,d,e,f}; do
	echo "$classe";
	echo $classe >> $fileout;
	for istanza in `ls ./istanze/$classe`; do
		#echo "---- INIZIO --------------------------------------------------------------------------";
		#echo "./istanze/$classe/$istanza";

		#bound=`./singola_simulazione.sh ./istanze/$classe/$istanza | grep "relmipgap =" | head -n 1 | tail -n 1 | cut -f 3 -d'=' | tr -d ' '`;
		#./singola_simulazione_lp.sh ./istanze/$classe/$istanza
		./MEWCP_dsdp $1 ./istanze/$classe/$istanza $fileout
		#echo $bound;
		#echo $bound >> $fileout;
		#echo "---- FINE --------------------------------------------------------------------------";

	done
done
