# SIMPLE SCRIPT TO RUN SAMPLE
Qsbumu="0.5"
GaussQsbymu="0.5"
RegMass="0.3"

L="32"
N="512"
a=`awk -v N=${N} -v L=${L} 'BEGIN { print  ( L / N ) }'`

kMin="0.5"
kMax="3.0"

dkBin="0.5"
dkStep="0.5"

kRefMin="0.5"
kRefMax="3.0"

DIR="NEW_CHEESE"

mkdir -p ${DIR}

NRap="100"
AS="0"
OUT="0"

echo ${L}
echo ${N}
echo ${a}

IF="input_EXAMPLE"

SEED="12847593"

./Simulation.exe -N ${N} -a ${a} -NRap $NRap} -o ${DIR} -AS ${AS} -OUT ${OUT} -kMin ${kMin} -kMax ${kMax} -kRefMin ${kRefMin} -kRefMax ${kRefMax} -dkBin ${dkBin} -dkStep ${dkStep} -k -RegMass ${RegMass} -Qsbymu ${Qsbumu} -GaussQsbymu ${GaussQsbymu} -if ${IF} -SEED ${SEED}

