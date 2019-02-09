# Build docker and deploy to docker hub
# Currently using local docker hub on cutter

docker build -t lama --no-cache --network=host .

docker tag lama cutter:5000/neil_lama
