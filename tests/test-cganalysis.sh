#!/bin/bash

source $HOME/.mummi/config.mummi.sh

python $MUMMI_APP/mummi_ras/online/cg/cganalysis.py --simname structure_iter00_000000000015 \
                                                    --noninteractive \
                                                    --path $MUMMI_ROOT/sims-cg/structure_iter00_000000000015 \
                                                    --fstype "mummi" \
                                                    --fbio "mummi" \
                                                    --logpath $MUMMI_ROOT/log \
                                                    -c \
                                                    --siminputs "" \
                                                    --fcount 10 \
                                                    --backend NOSIM \
                                                    --nprocs 1
