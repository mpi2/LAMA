# Build docker and deploy to docker hub
# Currently using local docker hub on cutter

#build
docker build -t neil_lama --no-cache --network=host .

#tag
docker tag neil_lama:latest cutter:5000/neil_lama

#push
docker push cutter:5000/neil_lama
