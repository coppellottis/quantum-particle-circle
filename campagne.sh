#!/bin/bash
# Uso:  ./campagne.sh plan         -> stampa N_t, NR, C, estremo caldo (non lancia nulla)
#       ./campagne.sh metropolis   -> campagna Metropolis (fase 2)
#       ./campagne.sh pt           -> campagna Parallel Tempering (fase 3)
# Va lanciato dalla cartella dove si trova ./quantum.

BHAT=2                                   # \hat\beta del bersaglio
NTS=(120 150 180 210 240 270 300 330 360 390)
PACC=0.5                                 # accettazione degli scambi voluta
NR_MAX=15                                # tetto sul numero di repliche
MEAS=1e6                                 # campionamenti per ogni simulazione
BHAT_MAX=${BHAT_MAX:-5}                  # estremo caldo: eta_max = BHAT_MAX/N_t  (si puo' cambiare: BHAT_MAX=4 ./campagne.sh plan)
THREADS=${THREADS:-$(nproc)}             # thread per la campagna Metropolis

pt_params() {   # $1 = N_t  ->  "NR C eta0 eta_top NR_grezzo eta_max"
python3 - "$1" "$BHAT" "$PACC" "$NR_MAX" "$BHAT_MAX" << 'PY'
import sys, math
from scipy.special import erfcinv
Nt = int(sys.argv[1]); bhat = float(sys.argv[2]); P = float(sys.argv[3]); nrmax = int(sys.argv[4])
etamax = float(sys.argv[5]) / Nt
eta0 = bhat / Nt
x = erfcinv(P)
C = 1 + 2 * math.sqrt(2) * x / math.sqrt(Nt)
raw = math.sqrt(Nt) / (2 * math.sqrt(2) * x) * math.log(etamax / eta0) + 1
NR = min(nrmax, max(2, math.ceil(raw)))
print("%d %.16g %.16g %.6g %.2f %.6g" % (NR, C, eta0, eta0 * C ** (NR - 1), raw, etamax))
PY
}

plan() {
  printf "%5s %4s %12s %11s %11s %9s %9s\n" N_t NR C eta_top bhat_top NR_grezzo eta_max
  for Nt in "${NTS[@]}"; do
    read -r NR C ETA0 ETOP RAW EMAX <<< "$(pt_params $Nt)"
    printf "%5d %4d %12.8f %11.5f %11.3f %9s %9s\n" $Nt $NR $C $ETOP "$(awk "BEGIN{print $ETOP*$Nt}")" $RAW $EMAX
  done
}

metropolis() {
  T=$(( THREADS < ${#NTS[@]} ? THREADS : ${#NTS[@]} ))
  { printf "met_bhat%s\n2\n%s\n%d\n" "$BHAT" "$BHAT" "${#NTS[@]}"
    printf "%d\n" "${NTS[@]}"
    printf "%s\n%d\n" "$MEAS" "$T"; } | ./quantum
}

pt() {
  for Nt in "${NTS[@]}"; do
    read -r NR C ETA0 ETOP RAW EMAX <<< "$(pt_params $Nt)"
    echo ">>> PT  N_t=$Nt  NR=$NR  C=$C  eta0=$ETA0  eta_top=$ETOP (bhat_top=$(awk "BEGIN{print $ETOP*$Nt}"))"
    printf "pt_N%d_bhat%s\n3\n%d\n%d\n%s\n%s\n%s\n" $Nt $BHAT $Nt $NR $C $ETA0 $MEAS | ./quantum
  done
}

case "$1" in
  plan) plan ;;
  metropolis) metropolis ;;
  pt) pt ;;
  *) echo "Uso: ./campagne.sh plan | metropolis | pt"; exit 1 ;;
esac
