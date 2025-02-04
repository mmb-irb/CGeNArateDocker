
# CGeNArate docker web services

Check that mongodb/import.sh has exec permissions, if not:

    chmod +x mongodb/import.sh

docker swarm init

docker network create --driver overlay --attachable dbnet

docker-compose build

export $(grep -v '^#' .env | xargs)
docker stack deploy -c docker-compose.yml my_stack


## Credits

David Farré-Gil, Genís Bayarri, Adam Hospital.

## Copyright & licensing

This website has been developed by the [MMB group](https://mmb.irbbarcelona.org) at the [IRB Barcelona](https://irbbarcelona.org).

© 2025 **Institute for Research in Biomedicine**

Licensed under the [**Apache License 2.0**](LICENSE).