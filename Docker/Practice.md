# Docker in Action

## Set up the environment
```
mdkir hello-docker
cd hello-docker
touch app.js
```

## Inside the app.js
```
console.log("Hello Docker!")
```

## Run the application in the terminal
```
node app.js
```

## Deploy the application without the Docker
* Start with an operating system 
* Then we need to install node which is an execution environment for javascript code 
* Next we need to copy our application files `app.js` 
* Finally we need to run `node app.js` 

If we are working with a sophisticated application, we would end up with a **complex release document that has to be precisely followed**. This is where docker comes to the rescue. We can write these instructions inside a docker file and let docker package up our application. 

## Inside the Dockerfile, document the deployment process
```
FROM node:alpine 
# Build a base image, which can be found in the Docker Hub. 
# :tag, we decide which linux distribution is used to build the node image.

COPY . /app
# Copy all the files in the current directory into the app directory

WORKDIR /app
# Set up the working directory

CMD node app.js
# Execute the command in the working directory
```

## Build the Docker in the terminal
```
docker build -t hello-docker .
# -t: give our image a tag named "hello-docker"
# .: the current directory stores the Docker file
```

## View all the images in the terminal
```
docker images
# TAG is used to version the docker images
# Each image has its own IMAGE ID
# The docker image size is 180 MB, as our image contains both the apline linux node and app.js
```

## Run the image in the terminal
```
docker run hello-docker
# It doesn't matter which directory we are in because this image contains all the files for running our application 
```

## Push the image into the Docker Hub
```
docker tag hello-docker:latest XYZ/hello-docker:latest
docker push XYZ/hello-docker:latest

# push a new image to this repository
# docker tag local-image-name:tag-name repo-name:tag-name
# docker push repo-name:tag-name

# push a new tag to this repository
# docker push repo-name:tag-name
```