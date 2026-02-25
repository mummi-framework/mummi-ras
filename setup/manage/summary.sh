#!/usr/bin/env bash

DIR="sims-cg"

RED="\e[31m"
BOLDRED="\e[1;31m"
BLUE="\e[34m"
BOLDBLUE="\e[1;34m"
YELLOW="\e[33m"
BOLDYELLOW="\e[1;33m"
GREEN="\e[32m"
BOLDGREEN="\e[1;32m"
PURPLE="\e[1;35m"
ENDCOLOR="\e[0m"

show_failure=0
parsable=0

print_usage() {
  echo "Usage: $0 [-s] [-h]"
  echo "	-s : show detailed failed createsims"
  echo "	-p : show list of failed createsims directory"
}

while getopts 'sph' flag; do
  case "${flag}" in
    s) show_failure=1 ;;
    p) parsable=1 ;;
    *) print_usage
       exit 1 ;;
  esac
done

if [ $parsable -eq 1 ]; then
    show_failure=1
fi

job_createsim=$(grep -A 3 jobs_createsim workspace/wfmanager.chk)
job_cg=$(grep -A 3 jobs_cg workspace/wfmanager.chk)
job_aa=$(grep -A 3 jobs_aa workspace/wfmanager.chk)
job_bm=$(grep -A 3 jobs_backmapping workspace/wfmanager.chk)

total_createsim=$(find $DIR/structure_iter*/ -name createsims_* | wc -l 2>/dev/null)
if [ $parsable -eq 0 ]; then
    if [ "$total_createsim" -eq "0" ]; then
        echo "Error: no simulations have finished in $DIR".
        exit
    else
        echo "> Processing $total_createsim simulations..."
    fi
fi

fail_createsim_list=$(find $DIR/structure_iter*/ -name createsims_failure)
fail_createsim=$(echo $fail_createsim_list | wc -l)
success_createsim=$(($total_createsim - $fail_createsim))
rate_createsim=$(($fail_createsim*100 / $total_createsim))

if [ $parsable -eq 0 ]; then
    echo -e "> MuMMI run in $PWD"
    echo -e " > Jobs                 : running / queued"
    echo -e " > Createsim            : $(echo $job_createsim | cut -d' ' -f7) / $(echo $job_createsim | cut -d' ' -f5)"
    echo -e "    ${BOLDBLUE}Total               : $total_createsim${ENDCOLOR}"
    echo -e "    ${BOLDGREEN}Successes           : $success_createsim${ENDCOLOR}"
    echo -e "    ${BOLDRED}Failures            : $fail_createsim${ENDCOLOR}"
    echo -e "${BOLDYELLOW}    Failures rate (%)   : $rate_createsim${ENDCOLOR}"
fi

if [ $show_failure -eq 1 ]; then
	if [ $fail_createsim -gt 0 ]; then
		if [ $parsable -eq 0 ]; then
            echo -e "    Failed createsims:"
        fi
		for s in $fail_createsim_list; do
			folder=$(echo "$s" | cut -d'/' -f1,2)
			if [ $parsable -eq 0 ]; then
                echo -e "      > ${BOLDRED}$folder${ENDCOLOR}"
			    #cause=$(grep -n "Failed run" $folder/failed_*/createsims.log | cut -d " " -f2,3,4,5,6,7,8 --complement | cut -d':' -f 1,3,4)
			    grep -Hn "Failed run" $folder/failed_*/createsims.log | cut -d " " -f2,3,4,5,6,7,8 --complement | cut -d':' -f 1,3,4 | while read -r line ; do
				    echo -e "        $line"
			    done
            else
                echo "$(realpath $folder)"
            fi
		done
	fi
fi

if [ $parsable -eq 0 ]; then
    stop_cg_list=$(find $DIR/structure_iter*/ -name cg_stop)
    fail_cg_list=$(find $DIR/structure_iter*/ -name cg_failure)
    success_cg=$(find $DIR/structure_iter*/ -name cg_success | wc -l 2>/dev/null)
    fail_cg=$(find $DIR/structure_iter*/ -name cg_failure | wc -l)
    stop_cg=$(find $DIR/structure_iter*/ -name cg_stop | wc -l)
    total_cg=$(find $DIR/structure_iter*/ -name cg_analysis.log | wc -l 2>/dev/null)
    
    if [ "$total_cg" -eq "0" ]; then
        rate_cg=0
    else
        rate_cg=$(($fail_cg*100 / $total_cg))
    fi

    echo -e " > CG                   : $(echo $job_cg | cut -d' ' -f7) / $(echo $job_cg | cut -d' ' -f5)"
    echo -e "    ${BOLDBLUE}Total               : $total_cg${ENDCOLOR}"
    echo -e "    ${BOLDGREEN}Successes           : $success_cg${ENDCOLOR}"
    echo -e "    ${PURPLE}Stopped             : $stop_cg${ENDCOLOR}"
    echo -e "    ${BOLDRED}Failures            : $fail_cg${ENDCOLOR}"
    echo -e "${BOLDYELLOW}    Failures rate (%)   : $rate_cg${ENDCOLOR}"

    if [ $fail_cg -gt 0 ]; then
        echo -e "    Failed cg:"
        for s in $fail_cg_list; do
            folder=$(echo "$s" | cut -d'/' -f1,2)
            #cause=$(grep "Failed run" $folder/failed_*/cg_analysis.out | cut -d " " -f1,2,3,4,5,6,7,8 --complement)
            echo -e "      > ${BOLDRED}$folder${ENDCOLOR}" # => $cause" 
        done
    fi
 
#    if [ $stop_cg -gt 0 ]; then
#        echo -e "    Stopped cg:"
#        for s in $stop_cg_list; do
#            folder=$(echo "$s" | cut -d'/' -f1,2)
#            #cause=$(grep "Failed run" $folder/failed_*/cg_analysis.out | cut -d " " -f1,2,3,4,5,6,7,8 --complement)
#            echo -e "      > ${PURPLE}$folder${ENDCOLOR}" # => $cause" 
#        done
#    fi

    echo -e " > AA          : $(echo $job_aa | cut -d' ' -f7) / $(echo $job_aa | cut -d' ' -f5)"
    echo -e " > Backmapping : $(echo $job_bm | cut -d' ' -f7) / $(echo $job_bm | cut -d' ' -f5)"
fi

