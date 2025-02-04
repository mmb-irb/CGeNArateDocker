#!/bin/bash

# Wait for MongoDB to be ready before importing the data
until mongosh "mongodb://${MONGO_INITDB_ROOT_USERNAME}:${MONGO_INITDB_ROOT_PASSWORD}@${DB_HOST}:${MONGO_PORT}/${MONGO_INITDB_DATABASE}?authSource=admin" --eval "db.stats()" >/dev/null 2>&1; do 
  echo "Waiting for MongoDB..."
  sleep 2
done

echo "MongoDB is ready. Importing ${DB_COLLECTION_WFS} collection..."
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
exit 0  # Ensure the container exits successfully