#!/bin/bash

mpiexec -n $1 python3 $MUMMI_APP/mummi_ras/online/cg/cg_analysis_executor.py --simname "" \
                                                    --noninteractive \
                                                    --path "" \
                                                    --fstype "mummi" \
                                                    --fbio "mummi" \
                                                    --siminputs "" \
                                                    --frameProcessBatchSize 100 \
                                                    --simbin "" \
                                                    --fcount 0\
                                                    --useGromacs

