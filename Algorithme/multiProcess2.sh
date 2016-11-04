#!/bin/bash

NUM=0
QUEUE=""
MAX_NPROC=1 # default - Pourrait devenir un paramètre d'entree de multiprocess2
REPLACE_CMD=0 # no replacement by default
function queue {
    QUEUE="$QUEUE $1"
    NUM=$(($NUM+1))
}
 
function regeneratequeue {
    OLDREQUEUE=$QUEUE
    QUEUE=""
    NUM=0
    for PID in $OLDREQUEUE
    do
        if [ -d /proc/$PID  ] ; then
            QUEUE="$QUEUE $PID"
            NUM=$(($NUM+1))
        fi
    done
}
 
function checkqueue {
    OLDCHQUEUE=$QUEUE
    for PID in $OLDCHQUEUE
    do
        if [ ! -d /proc/$PID ] ; then
            regeneratequeue # at least one PID has finished
            break
        fi
    done
}

ACO="mouvant.py" # Nom du script à executer
gcc -shared -Wl,-soname,mouvant -o mouvant.so -fPIC mouvant.c &
wait

for fileDat in  *_para.dat;
do
	nomSemaine=${fileDat:0:-9};
	filePara=$nomSemaine"_para.dat"
	fileTran=$nomSemaine"_tran.dat"
	cmd="python3 $ACO $filePara $fileTran" 
	$cmd &
	# DEFINE COMMAND END

	PID=$!
	queue $PID

	while [ $NUM -ge $MAX_NPROC ]; do
		checkqueue
		sleep 0.4
	done

done

    





exit 0
