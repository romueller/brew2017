#!/bin/bash

########## CONFIGURATION ##########

# Adjust to your machine

# Binaries
SWARM_OLD="..." # path to binary of old swarm version (breaking via script)
SWARM="..." # path to binary of new swarm version
GEFAST="..." # path to GeFaST binary
PYTHON="python2.7" # python command / path to binary


# Adjust only if you know what you are doing

# Scripts
SWARM_BREAKER="scripts/swarm_breaker.py"

###################################

# Inputs
FASTA=$1
RESULTS=$2
THRESHOLD=$3
TAXA=$4
METHOD=$5
THREADS=$6

# Perform Clustering
case ${METHOD} in

	# SWARM (v1.2.3)
	swarm-old) 
		"${SWARM_OLD}" -d ${THRESHOLD} -t ${THREADS} < "${FASTA}" > "${RESULTS}.tmp" 2> /dev/null; \
		"${PYTHON}" "${SWARM_BREAKER}" -b "${SWARM_OLD}" -f "${FASTA}" -s "${RESULTS}.tmp" -d ${THRESHOLD} 2> /dev/null > "${RESULTS}"; \
		rm "${RESULTS}.tmp"
	;;

	# SWARM (current version)
	swarm) "${SWARM}" -d ${THRESHOLD} -t ${THREADS} < "${FASTA}" > "${RESULTS}" 2> /dev/null
	;;
	
	# GeFaST (edit-distance)
	gefast-e) "${GEFAST}" "${FASTA}" -t ${THRESHOLD} -so "${RESULTS}" --per-worker ${THREADS} \
		> /dev/null 2>&1
	;;
	
	# GeFaST (edit-distance, fastidious, t + 1)
	gefast-ef1) "${GEFAST}" "${FASTA}" -t ${THRESHOLD} -so "${RESULTS}" --per-worker ${THREADS} \
		-sf --swarm-fastidious-threshold $((THRESHOLD + 1)) --swarm-num-threads-per-check ${THREADS} \
		> /dev/null 2>&1
	;;
	
	# GeFaST (edit-distance, fastidious, 2 * t)
	gefast-e2f) "${GEFAST}" "${FASTA}" -t ${THRESHOLD} -so "${RESULTS}" --per-worker ${THREADS} \
		-sf --swarm-num-threads-per-check ${THREADS} \
		> /dev/null 2>&1
	;;

	# GeFaST (scoring function)
	gefast-s) "${GEFAST}" "${FASTA}" -t ${THRESHOLD} -so "${RESULTS}" --per-worker ${THREADS} \
		--use-score \
		> /dev/null 2>&1
	;;

	# GeFaST (scoring function, fastidious, t + 1)
	gefast-sf1) "${GEFAST}" "${FASTA}" -t ${THRESHOLD} -so "${RESULTS}" --per-worker ${THREADS} \
		--use-score -sf --swarm-fastidious-threshold $((THRESHOLD + 1)) --swarm-num-threads-per-check ${THREADS} \
		> /dev/null 2>&1
	;;

	# GeFaST (scoring function, fastidious, 2 * t)
	gefast-s2f) "${GEFAST}" "${FASTA}" -t ${THRESHOLD} -so "${RESULTS}" --per-worker ${THREADS} \
		--use-score -sf --swarm-num-threads-per-check ${THREADS} \
		> /dev/null 2>&1
	;;

esac

