TO=600
MO=512
RESULTS_DIR=results
mkdir -p $RESULTS_DIR


# runlim -r $TO -s $MO  python3 -u main.py graph01  2>&1 | tee $RESULTS_DIR/graph01.log
# runlim -r $TO -s $MO  python3 -u main.py graph02  2>&1 | tee $RESULTS_DIR/graph02.log
# runlim -r $TO -s $MO  python3 -u main.py graph08  2>&1 | tee $RESULTS_DIR/graph08.log
# runlim -r $TO -s $MO  python3 -u main.py graph09  2>&1 | tee $RESULTS_DIR/graph09.log
# runlim -r $TO -s $MO  python3 -u main.py graph14  2>&1 | tee $RESULTS_DIR/graph14.log
# runlim -r $TO -s $MO  python3 -u main.py scen01  2>&1 | tee $RESULTS_DIR/scen01.log
# runlim -r $TO -s $MO  python3 -u main.py scen02  2>&1 | tee $RESULTS_DIR/scen02.log
# runlim -r $TO -s $MO  python3 -u main.py scen03  2>&1 | tee $RESULTS_DIR/scen03.log
# runlim -r $TO -s $MO  python3 -u main.py scen04  2>&1 | tee $RESULTS_DIR/scen04.log
runlim -r $TO -s $MO  python3 -u main.py scen11  2>&1 | tee $RESULTS_DIR/scen11.log