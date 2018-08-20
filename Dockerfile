# Use an official Python runtime as a parent image
FROM ubuntu:18.04

# Enable neurodebian
#RUN wget -O- http://neuro.debian.net/lists/bionic.de-md.libre | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
#RUN apt-key adv --recv-keys --keyserver hkp://pool.sks-keyservers.net:80 0xA5D32F012649A5A9

RUN apt update 
ENV DEBIAN_FRONTEND=noninteractive 
RUN apt install -y tzdata
RUN apt install -y r-base 
RUN apt install -y python-pip
RUN apt install nano
RUN apt install wget

#RUN apt install neurodebian
#RUN apt install fsl-5.0-core

WORKDIR /lama

# Copy the current directory contents into the container at /app
ADD . /lama

# Make port 80 available to the world outside this container
EXPOSE 80

# Install any needed packages specified in requirements.txt
RUN pip install --trusted-host pypi.python.org -r requirements.txt
