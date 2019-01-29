#!/bin/bash
###################################################
# Info: Shell script to manage different Docker options
# Author: Manu Molina
# Date: 23/01/2019
# version: 2.0
###################################################
checkAlias()
{
    data=$1
    if [ -z "$data" ]; then
        return "0"
    else
        return "1"
    fi
}

# Params:   alias value
#           
#
getAliasFromInput() {
    alias=`$fundialog --stdout --clear --title "name" --inputbox $1 0 0`
    return $alias
}

checkAliasAndWork() {

    if [ "$1" = "0" ]; then
        dialog --infobox "No alias, no Docker \n " 0 0
    else
        case $3 in
            2)  dialog --infobox "Your alias: $2 \n 
                    Your shared folder: $current \n
                    Running container: \n `docker run -dit --name=$2 -v $current:/root/src/karyon/shared --rm cgenomics/karyonpip:1.0`" 0 0;;
            #3)  dialog --infobox "Starting container:-\n `docker exec -it $2 /bin/bash`" 0 0;;
            3)  dialog --infobox "You are now inside $2\n" 0 0
                docker exec -it $2 /bin/bash
                break;;
            4)  dialog --infobox "Stoping container: \n `docker stop $2`" 0 0;;
            5)  dialog --infobox "Removing container: \n `docker rm $2`" 0 0;;
        esac
    fi
}


alias=
menu=
result=
current="$(pwd)/shared"


#
# set infinite loop
#
while true
do

### display main menu ###
fundialog=${fundialog=dialog}
menuitem=`$fundialog --stdout --backtitle "Karyon Docker Menu" \
--title "[ M A I N - M E N U ]" \
--menu "You can use the UP/DOWN arrow keys, the first \n\
letter of the choice as a hot key, or the \n\
number keys 1-8 to choose an option.\n\
Choose the OPTION" 0 0 0 \
1 "Build image" \
2 "Create container" \
3 "Start container " \
4 "Stop container" \
5 "Remove container" \
6 "View running processes" \
7 "View Docker images" \
8 "Quit"`

# make decsion 
case $menuitem in
    1)  dialog --infobox "Building image: \n " 0 0
        docker build -t cgenomics/karyonpip:1.0 .
        break;;
    2)  alias=`$fundialog --stdout --clear --title "name" --inputbox "Type the alias of the container you want to create:" 0 0`
        checkAlias $alias
        result=$?
        checkAliasAndWork $result $alias "2";;
    3)  alias=`$fundialog --stdout --clear --title "name" --inputbox "Type the alias of the container you want to start:" 0 0`
        checkAlias $alias
        result=$?
        checkAliasAndWork $result $alias "3";;
    4)  alias=`$fundialog --stdout --clear --title "name" --inputbox "Type the alias of the container you want to stop:" 0 0`
        checkAlias $alias
        result=$?
        checkAliasAndWork $result $alias "4";;
    5)  alias=`$fundialog --stdout --clear --title "name" --inputbox "Type the alias of the container you want to remove:" 0 0`
        checkAlias $alias
        result=$?
        checkAliasAndWork $result $alias "5";;        
    6)  dialog --infobox "Docker container details: \n `docker ps`" 0 0;;
    7)  dialog --infobox "Docker image details: \n `docker images`" 0 0;;
    8)  echo "Bye!"; break;;

esac

read -r -p "Press space to continue..." key

done

