 #!/usr/bin/env bash
set +x
set -u

#~ --- Created:			Fri Nov 15 08:40:25 CET 2019
#~ --- Last Modified:	Tue Dec  3 08:49:24 CET 2019

script_version="0.4.1"

scriptname="${0##*/}"

Color_Off='\033[0m'

# Reconstruction couleurs :::

# En "[0" on trouve avec une couleur d'écriture à "blanc" :

# 1 bright white
# 2 dark grey
# 3 normal
# 4 underlined
# 5 normal
# 6 normal
# 7 bg as fg and fg inversed
# 8 fg as bg (hence nothing visible ...)
# 9 strikethrough

# 11 -> 29 rien de bien notable ...


# 30 black					bk
# 31 red					rd
# 32 green					gn
# 33 orange (yellow?)		oe
# 34 blue					be
# 35 purple					pe
# 36 cyan					cn
# 37 Not Totally White?		nt

# 38 & 39 rien de bien notable

# From 40 -> 47 (inc.) sets background color
# same as 30 -> 37 but + 10 !!

# 48 -> 89 rien de bien notable

# 90 grey					gy
# 91 salmon					sn
# 92 green
# 93 yellow					yw
# 94 blue
# 95 purple
# 96 light blue				lb
# 97 bright white			bw

# 98 & 99 rien de bien notable

# From 100 -> 107 (inc.) sets background color
# same as 90 -> 97 but + 10 !!

# Nothing on my terminal after 107 ... up to 1048


# With modifiers : -> nbus [normal,bold,underlined,strikethrough]
# 0 normal
# 1 bold
# 2 nothing peculiar
# 3 nothing peculiar
# 4 underline
# 5 nothing peculiar
# 6 nothing peculiar
# 7 add grey background except for colors in ranges 30->37 and 90->97 where it shows black color on colored background
# 8 fg color == background color ???
# 9 strikethrough
# 10 Nothing peculiar
# 11 Nothing peculiar
# 12 Nothing peculiar
# 

function get_color {
   income_string=${1}
   modifier=""
   for ((nb=0;nb<${#income_string};nb++)); do
       if [ ${nb} -eq 0 ]; then
          #~ echo "String to parse ... ${income_string:nb:1}"
          case ${income_string:nb:1} in
              n) modifier="0" ;;
              b) modifier="1" ;;
              u) modifier="4" ;;
              s) modifier="9" ;;
              *) echo "Unrecognized modifier character ..."; exit 1 ;;
          esac
       fi
       if [ ${nb} -eq 1 ]; then
          #~ echo "String to parse ... ${income_string:nb:2}"
          case ${income_string:nb:2} in
              bk) color="30" ;;
              rd) color="31" ;;
              gn) color="32" ;;
              be) color="34" ;;
              pe) color="35" ;;
              cn) color="36" ;;
              gy) color="90" ;;
              yw) color="93" ;;
              bw) color="97" ;;
              *) echo "Unrecognized color characters ... ${income_string:nb:2}"; exit 1 ;;
          esac
       fi

   done
   color_global_var="\033[${modifier};${color}m"
return 0
} # end function get_color

#~ The End of All Things (op. cit.)
