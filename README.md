# Karyon Docker
To-Do
A one paragraph description about the container.

## Getting Started (To-Do)

If you haven't done so already, you'll need to "Clone or Download" this project.
```shell
git clone https://github.com/Gabaldonlab/redundans.git
```

## Prerequisities
In order to run this container you'll need docker installed.
* [Linux](https://docs.docker.com/linux/started/)

In addition, to work with the Karyon shell script menu it's necessary to install "dialog command".
The dialog command (https://bash.cyberciti.biz/guide/Bash_display_dialog_boxes) allows you to display a variety of questions or display messages using dialog boxes from a shell script.

```shell
sudo apt-get update
sudo apt-get install dialog
```

## Usage (To-Do)

In order to start to work with Karyon Docker you can do it within the shell script menu (Beginner Docker user) or directly with Docker commands (Advanced Docker user).

To use the shell script menu type:
```shell
sh menu.sh
```

These are the basic commands to work with Karyon Docker.(???)

### Building a new image
Building your docker container will take a few minutes the first time through. 
If you change your script or Dockerfile and need to rebuild, the subsequent times will be faster, since Docker builds things in layers.
```shell
docker build -t cgenomics/karyonpip:1.0 .
```

### Testing a new image (???)
<!-- To test the Docker container, we need to run our container with the data inside the shared folder.  -->
```shell

```

### Running the image

Every time you run a Docker container, it uses the Karyon Docker image to remake it.
When the container is closed, everything in it goes away! 
Also, for every one of those containers you need to name it with an alias.

```shell
docker run -dit --name=your-alias -v /your/karyon/location:/root/src/karyon/shared --rm cgenomics/karyonpip:1.0`" 0 0;;
```

### Launch Karyon script

```shell
docker exec -w /your/karyon/location/ karyon python scripts/karyon.py -l shared/your-data-file-one shared/your-data-file-two... -d shared/output -T
```

### Interactive session in a Karyon container

To start an interactive session in a container first you need to run the image with an alias.

```shell
docker exec -it your-alias /bin/bash
```

### Stoping a container

```shell
docker stop your-alias
```

### Removing a container

```shell
docker rm your-alias
```

### Showing running containers

```shell
docker ps
```

### Showing images

```shell
docker images
```

## Container Parameters (To-Do)

List the different parameters available to your container

<!-- ```shell
docker run give.example.org/of/your/container:v0.2.1 parameters
```

One example per permutation 

```shell
docker run give.example.org/of/your/container:v0.2.1
```

Show how to get a shell started in your container too

```shell
docker run give.example.org/of/your/container:v0.2.1 bash
``` -->

<!-- #### Environment Variables

* `VARIABLE_ONE` - A Description
* `ANOTHER_VAR` - More Description
* `YOU_GET_THE_IDEA` - And another -->

## Shared volume

* `/your/karyon/location/` - File location

## Useful File Locations

* `/your/karyon/location/shared` - Folder to share data with the container

* `/your/karyon/location/src/dependencies` - Folder with programs used by Karyon script


## Built With (To-Do)

* List the software v0.1.3
* And the version numbers v2.0.0
* That are in this container v0.3.2

## Find Us (To-Do)

* [GitHub](https://github.com/Gabaldonlab)
<!-- * [Quay.io](https://quay.io/repository/your/docker-repository) -->

## Support (To-Do)

If you have any issues or doubts, please send an email to ....

<!-- ## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the 
[tags on this repository](https://github.com/your/repository/tags).  -->

## Authors (To-Do)

* **Miguel Ángel Naranjo Ortiz** - *Pipeline work* - [MANaranjo](https://github.com/MANaranjo)
* **Manuel Molina Marín** - *Docker work* - [manumolina](https://github.com/manumolina)

## License (To-Do)

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

<!-- ## Acknowledgments

* People you want to thank
* If you took a bunch of code from somewhere list it here -->

