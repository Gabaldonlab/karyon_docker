###################################################
# Info: Shell script to manage different Docker options
# Author: Manu Molina
# Date: 23/01/2019
# version: 2.0
###################################################
checkalias()
{
	data=$1
	if [ -z "$data" ]; then
		return "0"
	else
		return "1"
	fi
}
alias=
menu=
result=
current="$(pwd)/shared"

echo "Select your option > "
echo "1º Build image"
echo "2º Create container"
echo "3º Start container "
echo "4º Stop container"
echo "5º Remove container"
echo "6º View running processes"
echo "7º View Docker images"
echo "8º Quit"
read menu

if [ "$menu" = "1" ]; then
	docker build -t cgenomics/karyonpip:1.0 .
elif [ "$menu" = "2" -o "$menu" = "3" -o "$menu" = "4" -o "$menu" = "5" ]; then
	if [ "$menu" = "2" ]; then
		echo "Type the alias of the container you want to create > "
	elif [ "$menu" = "3" ]; then
		echo "Type the alias of the container you want to start > "
	elif [ "$menu" = "4" ]; then
		echo "Type the alias of the container you want to stop > "
	fi
	
	read alias
	checkalias $alias
	result=$?
	if [ "$result" = "0" ]; then
		echo "No alias, no Docker"
	else
		echo "Your alias: $alias"
		echo "Your shared folder: $current"
		if [ "$menu" = "2" ]; then
			echo "Running container"
			docker run -dit --name=$alias -v $current:/root/src/karyon/shared --rm cgenomics/karyonpip:1.0
		elif [ "$menu" = "3" ]; then
			echo "Starting container"
			docker exec -it $alias /bin/bash
		elif [ "$menu" = "4" ]; then
			echo "Stoping container"
			docker stop $alias
		elif [ "$menu" = "5" ]; then
			echo "Removing container"
			docker rm $alias			
		fi
	fi
elif [ "$menu" = "6" ]; then
	docker ps
elif [ "$menu" = "7" ]; then
	docker images
elif [ "$menu" = "8" ]; then
	echo "Bye!"
else
	echo "No option, no Docker"
fi