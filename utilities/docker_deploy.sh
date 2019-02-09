# Build docker and deploy to docker hub
# Currently using local docker hub on cutter

#build
docker build -t lama_test --no-cache --network=host .

#tag
docker tag lama_test:latest cutter:5000/neil_lama

#push
docker push cutter:5000/neil_lama
