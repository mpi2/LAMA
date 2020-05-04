# Build docker and deploy to docker hub
# Currently using local docker hub on cutter
# Accepts one argument, the version tag of the image
# Run from lama root  $ lama/utilities/docker_deploy.sh

version=$1

#build
docker build -t neil_lama --no-cache --network=host .

#tag
docker tag neil_lama cutter:5000/neil_lama:${version}

#push
docker push cutter:5000/neil_lama:${version}
