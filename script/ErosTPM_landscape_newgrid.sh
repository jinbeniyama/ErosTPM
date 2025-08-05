#!/usr/bin/bash
# Code for Eros TPM
# NOTE: 
#   - Color-term correction was already done.
#   - runtpm doesnt work when low thermal intertia is too low.
#     add after '-p 3600'. when TI < 30.


# Object (for output files)
OBJ=433

# Input arguments =============================================================
## shape model
fOBJ=$1
## spin file
fSPIN=$2
## obs file
fOBS=$3
## ephem file
fEPH=$4
## output directory
OUTDIR=$5
mkdir -p ${OUTDIR}


# Fixed parameters ============================================================
## Emissivity
EPS=0.9
## Bond albedo from Domingue+2002
BondA=0.12
# Fixed parameters ============================================================


# Free parameters =============================================================
# Thermal intertia
TIarr=(10 20 30 40 50 60 70 80 90 100 
       120 140 160 180 200 220 240 260 280 300
       330 360 390 420 450 480 510 540 570 600)
 
# Roughness (crater) parameters
# We do not fix this. We use roughness in the literature for validation purposes.
CAarr=(90   90  90  90  90  90  90  90  90  90  90)
CRarr=(0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)
# Free parameters =============================================================


# Calculate number of combinations
echo N_TIarr = ${#TIarr[@]}
echo N_CAarr = ${#CAarr[@]}
echo N_CRarr = ${#CRarr[@]}


# ROOP for TI
for idxT in "${!TIarr[@]}"; 
do
    # ROOP for crater (roughness)
    for idxC in "${!CAarr[@]}"; 
    do
    TI="${TIarr[idxT]}" 
    CA="${CAarr[idxC]}" 
    CR="${CRarr[idxC]}"
        printf 'Inputs: (%s,%s,%s,%s)\n' "$fOBJ" "$fEPH" "$fSPIN" "$fOBS"
        printf 'Params: (%s,%s,%s,%s)\n' "$BondA" "$TI" "$CA" "$CR"

    # Number of points per rotation phase 
    TI_th=30
    if [ ${TI} -le ${TI_th} ]; then
        NPHASE=3600
    else
        NPHASE=360
    fi

    # Run TPM
    # wo/-f: Lagerros approximation 
    # w/-f : Full heat diffusion within the craters, much slower, but safer, more physical
    # Just output the command
    echo 'Command: echo'  ${fOBJ} ${fEPH} ${EPS} ${TI} ${BondA} ${CA} ${CR} '| runtpm -S' ${fSPIN} '-o' ${fOBS}  '-p' ${NPHASE} '>' ${OUTDIR}/tpmout_${OBJ}_brute_ti${TI}_ca${CA}_cr${CR}.dat '&'
    echo ""
    # Start TPM!
        echo ${fOBJ} ${fEPH} ${EPS} ${TI} ${BondA} ${CA} ${CR} | runlstpm -S ${fSPIN} -o ${fOBS} -p ${NPHASE} > ${OUTDIR}/lstpmout_${OBJ}_brute_ti${TI}_ca${CA}_cr${CR}.dat &
    done
	# Wait until all 10 processes are finished
	wait
done
