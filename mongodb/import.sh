#!/bin/bash

# Creates ssh keys for SGE and copy them into a shared volume
echo "Generating ssh keys for SGE..."
apt-get update && apt-get install openssh-client -y
mkdir /.ssh
ssh-keygen -t rsa -b 4096 -f /.ssh/id_rsa -N '' -q
cp /.ssh/id_rsa /.ssh/id_rsa.pub /keys
echo "Generating seed folder in /logs..."
mkdir -p /logs/seed

# Wait for MongoDB to be ready before importing the data
until mongosh "mongodb://${MONGO_INITDB_ROOT_USERNAME}:${MONGO_INITDB_ROOT_PASSWORD}@${DB_HOST}:${MONGO_PORT}/${MONGO_INITDB_DATABASE}?authSource=admin" --eval "db.stats()" >/dev/null 2>&1; do 
  echo "Waiting for MongoDB..."
  sleep 2
done

echo "MongoDB is ready. Dropping collection ${DB_COLLECTION_WFS}..."
mongosh --host ${DB_HOST} \
        --port ${MONGO_PORT} \
        -u ${MONGO_INITDB_ROOT_USERNAME} \
        -p ${MONGO_INITDB_ROOT_PASSWORD} \
        --authenticationDatabase admin \
        ${MONGO_INITDB_DATABASE} \
        --eval "db.getCollection('${DB_COLLECTION_WFS}').drop()"

echo "Importing ${DB_COLLECTION_WFS} collection..."
mongoimport --host ${DB_HOST} \
            --db ${MONGO_INITDB_DATABASE} \
            --port ${MONGO_PORT} \
            --username ${MONGO_INITDB_ROOT_USERNAME} \
            --password ${MONGO_INITDB_ROOT_PASSWORD} \
            --authenticationDatabase admin \
            --collection ${DB_COLLECTION_WFS} \
            --jsonArray \
            --file "/${DB_COLLECTION_WFS}.json"

echo "Importing ${DB_COLLECTION_PRJ} collection..."
mongoimport --host ${DB_HOST} \
            --db ${MONGO_INITDB_DATABASE} \
            --port ${MONGO_PORT} \
            --username ${MONGO_INITDB_ROOT_USERNAME} \
            --password ${MONGO_INITDB_ROOT_PASSWORD} \
            --authenticationDatabase admin \
            --collection ${DB_COLLECTION_PRJ} \
            --jsonArray \
            --file "/${DB_COLLECTION_PRJ}.json"

echo "Import completed successfully."

# Create a file to indicate that the seed is complete
touch /logs/seed/seed-complete

exit 0  # Ensure the container exits successfully