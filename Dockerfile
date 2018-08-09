# Use an official Python runtime as a parent image
FROM ubuntu:18.04

RUN apt update 
ENV DEBIAN_FRONTEND=noninteractive 
RUN apt install -y tzdata
#RUN apt install tzdata -y 
RUN apt install -y r-base 
RUN  apt install -y python-pip

WORKDIR /lama

# Copy the current directory contents into the container at /app
ADD . /lama

# Make port 80 available to the world outside this container
EXPOSE 80

# Install any needed packages specified in requirements.txt
RUN pip install --trusted-host pypi.python.org -r requirements.txt
