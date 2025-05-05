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
TIarr=(0 10 20 30 40 50 60 70 80 90 100 
       150 200 250 300 350 400 450 500 550 600 
       650 700 750 800 850 900 950 1000
       1100 1200 1300 1400 1500 1600 1700 1800 1900 2000
       2100 2200 2300 2400 2500)

# Roughness (crater) parameters
# We do not fix this. We use roughness in the literature for validation purposes.
# The same as Hung+2022, PSJ
CAarr=(0    30  40  41  50  60  70  88  89  90)
CRarr=(0.0 0.3 0.7 0.9 0.5 0.9 0.7 0.5 0.7 0.9)
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
    echo 'Command: echo'  ${fOBJ} ${fEPH} ${EPS} ${TI} ${BondA} ${CA} ${CR} '| runtpm -f -S' ${fSPIN} '-o' ${fOBS}  '-p' ${NPHASE} '>' ${OUTDIR}/tpmout_${OBJ}_brute_ti${TI}_ca${CA}_cr${CR}.dat '&'
    echo ""
    # Start TPM!
        echo ${fOBJ} ${fEPH} ${EPS} ${TI} ${BondA} ${CA} ${CR} | runtpm -f -S ${fSPIN} -o ${fOBS} -p ${NPHASE} > ${OUTDIR}/tpmout_${OBJ}_brute_ti${TI}_ca${CA}_cr${CR}.dat &
    done
	# Wait until all 10 processes are finished
	wait
done
