#!/bin/bash

# This script takes H-representation of network after mplrs redund and iterates with mplrs project over it.

# Help page

usage="$(basename "$0") [-h] [-i FILE] [-n INT] [-c INT] [-m STR] -- This script takes H-representation of network after mplrs redund and iterates with mplrs project over it. Saves output in directory where executed.

where:
    -h  show this help text
    -i  inputfile
    -n  number of iterations
    -c  number of used cores for mplrs (default: 3)
    -m  path to mplrs (default: mplrs --> will search on its own)"

while getopts ':hi:n:c:m:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    i) inputfile=$OPTARG
       ;;
    n) number=$OPTARG
       ;;
    c) n_cores=$OPTARG
       ;;
    m) path_to_mplrs=$OPTARG
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))

start_time=$SECONDS

# mplrs project function
#mplrs_project () {
#    mpirun -np $n_cores $path_to_mplrs $inputfile.tmp -maxbuf 500 > $inputfile.projected;
#    mv -f $inputfile.projected $inputfile.tmp;
#}

# n_cores default is 3
if [ -z "$n_cores" ]; then n_cores=3; fi
echo "Running mplrs with $n_cores cores."

# path_to_mplrs default is mplrs
if [ -z "$path_to_mplrs" ]; then path_to_mplrs="mplrs"; fi

# prepare inputfile
cp $inputfile $inputfile.tmp

# loop over tmp file
for ((n=0;n<$number;n++)); do
    #echo "Run $(($n + 1))";
    #mplrs_project ()
    mpirun -np $n_cores $path_to_mplrs $inputfile.tmp -maxbuf 500 > $inputfile.projected;
    mv -f $inputfile.projected $inputfile.tmp;
done

# rename finished file again
mv -f $inputfile.tmp $inputfile.projected


elapsed=$(( SECONDS - start_time ))
echo "Ran in $elapsed seconds."



#https://phoenixnap.com/kb/bash-function




#n = len(smatrix[0]) - len(ex_reactions) + 1 # +1 because last one is the conversion from H to V-Representation

# if we want the last H-representation
    #n = len(smatrix[0]) - len(ex_reactions)

#def mplrs_project_3(smatrix, ex_reactions, n_cores, path_mplrs, core_name, n):
#    """
#    Performs mplrs project.
#    """
#    if n == 0:
#        if os.path.isfile(core_name + "_postredund.ine"):
#           os.replace(core_name + "_postredund.ine", core_name + ".projected")
#        return "Finished"
#    # mplrs project
#    else:
#        cmd = ["mpirun", "-np", str(n_cores), path_mplrs, core_name + "_postredund.ine", core_name + ".projected", "-maxbuf", str(500)]
#        subprocess.run(cmd, capture_output=True)
#        n -= 1
#        if n > 0:
#            cmd = ["mpirun", "-np", str(n_cores), path_mplrs, core_name + ".projected", core_name + "_postredund.ine", "-maxbuf", str(500)]
#            subprocess.run(cmd, capture_output=True)
#            n -= 1
#    mplrs_project_3(smatrix, ex_reactions, n_cores, path_mplrs, core_name, n)
