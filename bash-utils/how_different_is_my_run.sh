#!/bin/bash
set +x
set -u

script_version="0.0.5"

scriptname="${0##*/}"

Help()
{
   # Display Help
   echo -e "${IPurple}This function requests two arguments${Color_Off}"
   echo
   echo "Syntax: ${scriptname} [path_to_sources] [path_to_root_iLOVECLIM]"
   echo
}

ErrorNotDirectory()
{
   # Display Error
   echo -e "   ---   ${IRed}SCRIPT ERROR${Color_Off}   ---   "
   echo -e "The two arguments ${IRed}must be${Color_Off} directories"
   echo
   echo "[path_to_sources]        ... should contain a working version of the sources to compare against a reference -|"
   echo "[path_to_root_iLOVECLIM] ... should contain the root level of an iLOVECLIM installation                     -|"
   echo
}

# ====================   COLORIZATION    ==================== #

. Colors_Module.sh

color_global_var=""

get_color "bgn"
IGreen=${color_global_var}
get_color "ngy"
Gray=${color_global_var}
get_color "bpe"
IPurple=${color_global_var}
get_color "byw"
IYellow=${color_global_var}
get_color "brd"
IRed=${color_global_var}

blanks="                                                               "

echo
echo

if [ $# -ne 2 ]; then
  Help
  exit 1
fi

if [ -d "${1}" ] && [ -d "${2}" ] ; then
    pathanalysis=$(realpath ${1})
    pathrefcode=$(realpath ${2})
    echo "Starting analysis with ${pathanalysis} against the subdirectories of ${pathrefcode}";
else
    ErrorNotDirectory
    exit 1
fi

n=30

for fich in $(ls ${pathanalysis}/*.[fFh]*); do
    namefich=$(basename ${fich})
    reffich=$(find ${pathrefcode} -type f -regextype posix-egrep -regex ".*/(ecbilt|ecbilt_clio|clio|vecode|lbm|ludus-code|ludus_base_libs|ocycc|veccarb|medusa2|IO_NC|new_iceberg|downscaling|caraib_cpl|dgc-waffle|isoatm|lcm2ism)/.+" -name "${namefich}" | grep -v tryout_f90)
    if [ -z ${reffich} ]; then
        echo -e "\r Could not find reference location for : ${IRed}${namefich}${Color_Off} | New file ? \t \t \t"
    else
        diff -Bbq ${fich} ${reffich} > /dev/null
        res=$?
        if [ ${res} -eq 0 ]; then
           echo -en "\r ${Gray}${namefich} is identical to reference${Color_Off} ${blanks}"
        else
           echo -en "\r ${blanks} ${blanks}"
           refloc=${reffich#${pathrefcode}}
           refloc=${refloc%${namefich}}
           refloc=${refloc#/}
           refloc=${refloc%/}
           nmf=$namefich$(printf '%*s' "$n" "")                     # pad with `n` spaces.
           nmf=$(echo "${nmf}"|grep -Eo "^.{1,$n}")            # limit length to `n`
           echo -e "\r ${IPurple}${nmf}${Color_Off} is different from reference in ${IYellow}${refloc}${Color_Off}"
           # cp -vf ${fich} ${reffich}
        fi
    fi
done 

echo -en "\r ${blanks} ${blanks}"
echo 
